#' Construct a functional basis for the modulating component W
#'
#' Define the finite-dimensional space of functions \eqn{W: \Theta \to
#' \mathbb{R}^{p \times d}} from which the modulating component of the
#' AMM canonical form is drawn. The basis is evaluated at the anchor
#' \eqn{\theta_0} and at the working values of the population reference
#' \eqn{\theta_{\mathrm{ref}}} during fitting.
#'
#' For \eqn{p > 1}, the package-provided basis types (\code{"polynomial"}
#' and \code{"bspline"}) implement the \emph{separable} case: each
#' coordinate of \eqn{\theta_{\mathrm{ref}}} contributes an independent
#' block of basis functions, concatenated in dimension order. Non-separable
#' (cross-coupling) bases are planned for a future release together with
#' the Gaussian-process extension; the current API is forward-compatible.
#'
#' @param type Character scalar identifying the type of basis. One of
#'   \code{"polynomial"}, \code{"bspline"}, \code{"user"}.
#' @param degree Positive integer with the polynomial degree
#'   (\code{type = "polynomial"}) or the B-spline degree
#'   (\code{type = "bspline"}). Ignored for \code{type = "user"}.
#'   When \code{NULL} (default), resolves to \code{1L} for
#'   \code{type = "polynomial"} and \code{3L} for \code{type = "bspline"}
#'   (the conventional cubic-spline default).
#' @param knots Numeric vector of interior knots used by
#'   \code{splines::bs} when \code{type = "bspline"}. Required for
#'   \code{type = "bspline"} when \code{df} is not supplied.
#' @param df Positive integer with the number of basis functions used by
#'   \code{splines::bs} when \code{type = "bspline"}. Required when
#'   \code{knots} is not supplied. Mutually exclusive with \code{knots}.
#' @param boundary_knots Numeric vector of length two giving the
#'   univariate boundary knots of the B-spline expansion. Required for
#'   \code{type = "bspline"} when the basis will be wired through the
#'   Stan-side fit path (Path 1, sub-phase 8.3.8): the Cox--de Boor
#'   recursion in Stan needs a fixed boundary, while \code{splines::bs}
#'   would otherwise infer them from each evaluation point and produce
#'   knot vectors that drift between Hamiltonian Monte Carlo steps.
#'   Must satisfy \code{boundary_knots[1] < min(knots)} and
#'   \code{boundary_knots[2] > max(knots)} when \code{knots} is given.
#'   Ignored for \code{type != "bspline"}.
#' @param basis_fn A function that takes a numeric vector of length
#'   \code{p} (the value of theta_ref) and returns a numeric vector of
#'   length equal to the dimension of the basis. Required for
#'   \code{type = "user"}; ignored otherwise. The function must be pure
#'   (no side effects) and must return the same dimension for any input
#'   of length \code{p}.
#' @param dim Positive integer giving the dimension of the basis.
#'   Required for \code{type = "user"}; computed automatically for
#'   \code{type = "polynomial"} and \code{type = "bspline"}.
#' @param p Optional positive integer giving the dimension of
#'   \code{theta_ref}. When supplied, the basis is eagerly materialized
#'   at construction time: the package probes the evaluator, records
#'   \code{dim}, and computes per-coordinate block indices
#'   (\code{block_indices}) for separable types. When \code{NULL}
#'   (default), materialization is deferred to
#'   \code{materialize_W_basis()}.
#'
#' @return An object of class \code{W_basis} with components
#'   \code{type}, \code{degree}, \code{knots}, \code{df},
#'   \code{boundary_knots}, \code{dim}, \code{evaluator}, \code{p}, and
#'   \code{block_indices}.
#'   \code{evaluator} is a function that maps a numeric value of
#'   theta_ref (length \code{p}) to a numeric vector of length \code{dim}.
#'   \code{block_indices}, when populated, is a list of length \code{p}
#'   whose \eqn{k}-th entry contains the row indices of the basis output
#'   corresponding to coordinate \eqn{k} of theta_ref; it is \code{NULL}
#'   for \code{type = "user"} (separability of user bases cannot be
#'   inferred automatically).
#'
#' @details
#' Unlike the additive component \eqn{a(x)} and the multiplicative
#' component \eqn{b(x)}, which are defined on the covariate space
#' \eqn{\mathcal{X}} and are declared via standard R formulas, the
#' modulating component \eqn{W(\theta)} is defined on the parameter
#' space \eqn{\Theta}. The user therefore declares its basis functionally
#' rather than via \code{\link[stats]{model.matrix}}: the basis describes
#' how W depends on the value of theta_ref, not on the observed
#' covariates.
#'
#' Given a value of theta_ref of length \eqn{p}, the basis evaluator
#' returns a numeric vector of length \code{dim} encoding the basis
#' functions at that point. The Stan code generator combines these
#' values with the basis coefficients (one matrix of size
#' \eqn{p \times d} per basis function) to assemble \eqn{W(\theta)} at
#' every Hamiltonian Monte Carlo step.
#'
#' For \code{type = "polynomial"}, the basis includes monomials of
#' degrees 1 through \code{degree} of every coordinate of theta_ref,
#' arranged \strong{block-by-coordinate}: for \eqn{p > 1} the output is
#' \eqn{(\theta_1, \theta_1^2, \ldots, \theta_1^{\mathrm{degree}},
#' \theta_2, \theta_2^2, \ldots, \theta_p^{\mathrm{degree}})}.
#' Cross-terms between coordinates are excluded by default to keep the
#' basis dimension manageable; they can be added by supplying a
#' user-defined basis with \code{type = "user"}.
#'
#' For \code{type = "bspline"}, the basis is a B-spline expansion of a
#' single coordinate of theta_ref via \code{splines::bs}. When theta_ref
#' has more than one coordinate, the B-spline basis is applied to each
#' coordinate independently and \strong{concatenated block-by-coordinate},
#' matching the polynomial convention.
#'
#' @section Methodological notes:
#' The anchor point \eqn{\theta_0} at which \eqn{W(\theta_0) = 0} is a
#' parametrization device, not an inferential statement about the
#' posterior of theta_ref. Choosing a different anchor changes the
#' parametrization but not the data-generating model, provided the basis
#' spans the same function space. See assumption (C4) of Block 1 and the
#' anchoring discussion in Section 6.7 (Theorem 1E).
#'
#' Identifiability of \eqn{W} as a function on \eqn{\Theta} requires the
#' prior on theta_ref to assign positive measure to a connected open
#' subset of \eqn{\Theta} (assumption (BAY-1) of Theorem 1E). The basis
#' declared here is the finite-dimensional ambient space; the prior
#' over the basis coefficients is configured via
#' \code{\link{gdpar_prior}}.
#'
#' @section Dependencies:
#' This function uses \code{\link[splines]{bs}} from the \pkg{splines}
#' package (a base R package) when \code{type = "bspline"}.
#'
#' @references
#' See \code{vignette("v01_amm_identifiability", package = "gdpar")},
#' Section 6.7 (Theorem 1E) for the identifiability of W as a function
#' on the prior support.
#'
#' @examples
#' wb_poly <- W_basis(type = "polynomial", degree = 2)
#' print(wb_poly)
#' wb_poly$evaluator(0.5)
#'
#' wb_poly_mv <- W_basis(type = "polynomial", degree = 2, p = 2L)
#' print(wb_poly_mv)
#' wb_poly_mv$block_indices
#'
#' wb_user <- W_basis(
#'   type     = "user",
#'   basis_fn = function(theta) c(theta, theta^2, sin(theta)),
#'   dim      = 3
#' )
#' wb_user$evaluator(0.5)
#'
#' @seealso \code{\link{amm_spec}}, \code{\link{gdpar_prior}},
#'   \code{\link{as_per_k}}
#' @export
W_basis <- function(type = c("polynomial", "bspline", "user"),
                    degree = NULL, knots = NULL, df = NULL,
                    boundary_knots = NULL,
                    basis_fn = NULL, dim = NULL, p = NULL) {
  type <- match.arg(type)
  if (!is.null(p)) {
    assert_count(p, "p")
    p <- as.integer(p)
  }
  if (type == "polynomial") {
    if (is.null(degree)) degree <- 1L
    assert_count(degree, "degree")
    degree_int <- as.integer(degree)
    evaluator <- function(theta) {
      theta_num <- as.numeric(theta)
      p_theta <- length(theta_num)
      if (p_theta == 1L) {
        return(theta_num ^ seq_len(degree_int))
      }
      unlist(lapply(seq_len(p_theta),
                    function(k) theta_num[k] ^ seq_len(degree_int)))
    }
    obj <- list(type = type, degree = degree_int,
                knots = NULL, df = NULL,
                boundary_knots = NULL,
                dim = NA_integer_,
                evaluator = evaluator,
                p = NULL,
                block_indices = NULL)
  } else if (type == "bspline") {
    if (is.null(degree)) degree <- 3L
    if (is.null(knots) && is.null(df)) {
      gdpar_abort(
        paste(
          "For type = 'bspline', either 'knots' (interior knot vector)",
          "or 'df' (degrees of freedom) must be supplied."
        ),
        class = "gdpar_input_error"
      )
    }
    if (!is.null(knots) && !is.null(df)) {
      gdpar_abort(
        "Arguments 'knots' and 'df' are mutually exclusive.",
        class = "gdpar_input_error"
      )
    }
    assert_count(degree, "degree")
    if (!is.null(boundary_knots)) {
      if (!is.numeric(boundary_knots) || length(boundary_knots) != 2L ||
          any(!is.finite(boundary_knots)) ||
          boundary_knots[1L] >= boundary_knots[2L]) {
        gdpar_abort(
          paste(
            "boundary_knots must be a numeric vector of length 2 with",
            "boundary_knots[1] < boundary_knots[2] and both finite."
          ),
          class = "gdpar_input_error",
          data = list(received = boundary_knots)
        )
      }
      if (!is.null(knots) &&
          (min(knots) <= boundary_knots[1L] ||
           max(knots) >= boundary_knots[2L])) {
        gdpar_abort(
          paste(
            "Interior 'knots' must lie strictly inside 'boundary_knots';",
            "received knots range [", format(min(knots)),
            ", ", format(max(knots)), "] vs boundary [",
            format(boundary_knots[1L]), ", ",
            format(boundary_knots[2L]), "]."
          ),
          class = "gdpar_input_error",
          data = list(knots_range = range(knots),
                      boundary_knots = boundary_knots)
        )
      }
    }
    require_suggested("splines", "evaluate B-spline bases")
    .knots_local <- knots
    .df_local <- df
    .boundary_local <- boundary_knots
    .degree_local <- as.integer(degree)
    evaluator <- function(theta) {
      basis_mats <- lapply(theta, function(tk) {
        b <- if (!is.null(.knots_local)) {
          if (!is.null(.boundary_local)) {
            splines::bs(tk, knots = .knots_local, degree = .degree_local,
                        Boundary.knots = .boundary_local)
          } else {
            splines::bs(tk, knots = .knots_local, degree = .degree_local,
                        Boundary.knots = range(c(.knots_local, tk)))
          }
        } else {
          if (!is.null(.boundary_local)) {
            splines::bs(tk, df = .df_local, degree = .degree_local,
                        Boundary.knots = .boundary_local)
          } else {
            splines::bs(tk, df = .df_local, degree = .degree_local)
          }
        }
        as.numeric(b)
      })
      unlist(basis_mats)
    }
    obj <- list(type = type, degree = as.integer(degree),
                knots = knots, df = df,
                boundary_knots = boundary_knots,
                dim = NA_integer_,
                evaluator = evaluator,
                p = NULL,
                block_indices = NULL)
  } else {
    if (!is.function(basis_fn)) {
      gdpar_abort(
        "For type = 'user', 'basis_fn' must be a function.",
        class = "gdpar_input_error"
      )
    }
    assert_count(dim, "dim")
    evaluator <- basis_fn
    obj <- list(type = type, degree = NA_integer_,
                knots = NULL, df = NULL,
                boundary_knots = NULL,
                dim = as.integer(dim),
                evaluator = evaluator,
                p = NULL,
                block_indices = NULL)
  }
  class(obj) <- c("W_basis", "list")
  if (!is.null(p)) {
    obj <- materialize_W_basis(obj, p)
  }
  obj
}

#' Internal: finalize a W_basis at a given theta_ref dimension
#'
#' Records the total basis dimension \code{dim}, the theta_ref dimension
#' \code{p}, and (for separable types) the per-coordinate block indices
#' \code{block_indices}. For polynomial and B-spline bases the dimension
#' is derived analytically from the constructor arguments
#' (\code{degree}, \code{df}, \code{knots}), avoiding any runtime probing
#' of the evaluator; this both is more robust (no boundary-knot
#' degeneracy warnings from \code{splines::bs} at a single probe point)
#' and more efficient. For \code{type = "user"} the evaluator is probed
#' at \code{rep(0.5, p)} and the user-declared \code{dim} is validated
#' against the actual output length; \code{block_indices} is left
#' \code{NULL} because separability of an arbitrary user-supplied
#' evaluator cannot be inferred automatically.
#'
#' @param wb An object of class \code{W_basis}.
#' @param p Length of theta_ref.
#' @return A modified \code{W_basis} object with \code{dim}, \code{p},
#'   and (when computable) \code{block_indices} populated.
#' @keywords internal
#' @noRd
materialize_W_basis <- function(wb, p) {
  assert_inherits(wb, "W_basis", "wb")
  assert_count(p, "p")
  p <- as.integer(p)
  if (!is.null(wb$p) && !is.na(wb$p) && wb$p != p) {
    gdpar_abort(
      sprintf(
        paste(
          "W_basis was constructed with p = %d but materialize_W_basis",
          "was called with p = %d. These must match."
        ),
        wb$p, p
      ),
      class = "gdpar_input_error",
      data = list(wb_p = wb$p, requested_p = p)
    )
  }
  if (wb$type == "polynomial") {
    per_k_dim <- as.integer(wb$degree)
    wb$dim <- per_k_dim * p
    wb$block_indices <- lapply(
      seq_len(p),
      function(k) ((k - 1L) * per_k_dim + 1L):(k * per_k_dim)
    )
  } else if (wb$type == "bspline") {
    per_k_dim <- if (!is.null(wb$df)) {
      as.integer(wb$df)
    } else {
      as.integer(length(wb$knots) + wb$degree)
    }
    wb$dim <- per_k_dim * p
    wb$block_indices <- lapply(
      seq_len(p),
      function(k) ((k - 1L) * per_k_dim + 1L):(k * per_k_dim)
    )
  } else {
    probe <- rep(0.5, p)
    out <- tryCatch(
      wb$evaluator(probe),
      error = function(e) {
        gdpar_abort(
          sprintf(
            paste(
              "Failed to evaluate user W_basis at probe theta_ref of",
              "length %d: %s"
            ),
            p, conditionMessage(e)
          ),
          class = "gdpar_input_error",
          data = list(probe = probe)
        )
      }
    )
    if (!is.numeric(out)) {
      gdpar_abort(
        "The user W_basis evaluator must return a numeric vector.",
        class = "gdpar_input_error",
        data = list(returned_class = class(out))
      )
    }
    if (length(out) != wb$dim) {
      gdpar_abort(
        sprintf(
          paste(
            "User-supplied W_basis declares dim = %d but evaluator returned",
            "length %d for theta_ref of length %d."
          ),
          wb$dim, length(out), p
        ),
        class = "gdpar_input_error"
      )
    }
    wb$block_indices <- NULL
  }
  wb$p <- p
  wb
}

#' Split a separable multivariate W_basis into per-coordinate sub-bases
#'
#' For separable multivariate W bases (polynomial and bspline with
#' \eqn{p > 1}), return a list of length \code{p} whose \eqn{k}-th entry
#' is a univariate \code{W_basis} describing the contribution of
#' coordinate \eqn{k} of \code{theta_ref}. Useful for introspection,
#' per-coordinate diagnostics (e.g., the per-\eqn{k} CP/NCP toggle of
#' the pre-flight pipeline), and downstream methods that operate on one
#' coordinate at a time.
#'
#' @param wb A \code{W_basis} object with \code{p} populated (either
#'   constructed with \code{p = ...} or processed by
#'   \code{materialize_W_basis()}).
#'
#' @return A list of length \code{wb$p}. Each entry is itself a
#'   \code{W_basis} of the same \code{type} as \code{wb}, materialized
#'   with \code{p = 1L}. For \code{type = "user"} the function returns
#'   \code{NULL}: separability of a user-supplied evaluator cannot be
#'   inferred automatically. A warning is emitted in that case.
#'
#' @details
#' Polynomial and bspline bases are separable by construction: the
#' total basis is the concatenation, in coordinate order, of the
#' univariate bases applied independently to each coordinate of
#' \code{theta_ref}. \code{as_per_k()} simply reconstructs that
#' decomposition as explicit \code{W_basis} objects, each materialized
#' with \code{p = 1L}.
#'
#' This decomposition is the natural input to the per-coordinate
#' CP/NCP toggle of the pre-flight pipeline (Path B'), which operates
#' on each coordinate of \code{theta_ref} independently when \code{p > 1}.
#'
#' For \code{type = "user"} the separability is unverifiable and the
#' function returns \code{NULL}. Callers that need per-coordinate
#' decompositions for user bases must construct them explicitly.
#'
#' @examples
#' wb <- W_basis(type = "polynomial", degree = 2, p = 3L)
#' subs <- as_per_k(wb)
#' length(subs)
#' subs[[1]]$dim
#'
#' @seealso \code{\link{W_basis}}
#' @export
as_per_k <- function(wb) {
  assert_inherits(wb, "W_basis", "wb")
  if (is.null(wb$p) || is.na(wb$p)) {
    gdpar_abort(
      paste(
        "as_per_k() requires a materialized W_basis with 'p' populated.",
        "Either construct W_basis(..., p = ...) or call",
        "materialize_W_basis(wb, p)."
      ),
      class = "gdpar_input_error"
    )
  }
  if (wb$type == "user") {
    warning(
      "as_per_k() returns NULL for user-supplied W bases: ",
      "separability cannot be inferred automatically.",
      call. = FALSE
    )
    return(NULL)
  }
  if (wb$type == "polynomial") {
    out <- lapply(
      seq_len(wb$p),
      function(k) W_basis(type = "polynomial",
                          degree = wb$degree, p = 1L)
    )
    return(out)
  }
  if (wb$type == "bspline") {
    out <- lapply(seq_len(wb$p), function(k) {
      W_basis(type = "bspline",
              degree = wb$degree,
              knots = wb$knots,
              df = wb$df,
              boundary_knots = wb$boundary_knots,
              p = 1L)
    })
    return(out)
  }
  NULL
}

#' Print method for W_basis objects
#'
#' @param x An object of class \code{W_basis}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.W_basis <- function(x, ...) {
  cat("<W_basis>\n")
  cat("  type   : ", x$type, "\n", sep = "")
  if (!is.na(x$degree)) cat("  degree : ", x$degree, "\n", sep = "")
  if (!is.null(x$knots)) {
    cat("  knots  : ", paste(format(x$knots), collapse = ", "),
        "\n", sep = "")
  }
  if (!is.null(x$df)) cat("  df     : ", x$df, "\n", sep = "")
  if (!is.null(x$boundary_knots)) {
    cat("  boundary_knots: ",
        paste(format(x$boundary_knots), collapse = ", "),
        "\n", sep = "")
  }
  if (!is.na(x$dim)) cat("  dim    : ", x$dim, "\n", sep = "")
  if (!is.null(x$p) && !is.na(x$p)) {
    label <- if (x$p == 1L) {
      " (univariate)"
    } else if (!is.null(x$block_indices)) {
      " (multivariate, separable)"
    } else {
      " (multivariate)"
    }
    cat("  p      : ", x$p, label, "\n", sep = "")
  }
  if (!is.null(x$block_indices)) {
    cat("  block_indices: ",
        paste(vapply(x$block_indices,
                     function(idx) {
                       if (length(idx) == 1L) {
                         as.character(idx)
                       } else {
                         sprintf("%d-%d", min(idx), max(idx))
                       }
                     },
                     character(1L)),
              collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

#' Internal: resolve interior knots of a B-spline basis from df + boundary
#'
#' Sub-phase 8.3.8 (2026-05-22). When the user supplies a B-spline
#' \code{W_basis} via \code{df} (degrees of freedom) without explicit
#' \code{knots}, the Stan-side fit path still needs a fixed knot vector
#' to feed Cox--de Boor recursion (the basis must be deterministic
#' across Hamiltonian Monte Carlo steps). This helper derives the
#' interior knots by equally spacing them inside \code{boundary_knots},
#' yielding \code{df - degree} interior knots. Returns the input
#' \code{W_basis} with \code{knots} populated; passthrough when
#' \code{knots} is already set or when \code{type != "bspline"}.
#'
#' @keywords internal
#' @noRd
.gdpar_resolve_bspline_knots <- function(W) {
  if (W$type != "bspline") return(W)
  if (!is.null(W$knots)) return(W)
  if (is.null(W$df) || is.null(W$boundary_knots)) {
    gdpar_abort(
      paste(
        "Cannot resolve B-spline knots: when 'knots' is NULL the basis",
        "requires both 'df' and 'boundary_knots' to derive interior",
        "knots equally spaced inside the boundary."
      ),
      class = "gdpar_input_error",
      data = list(df = W$df, boundary_knots = W$boundary_knots)
    )
  }
  n_int <- as.integer(W$df) - as.integer(W$degree)
  if (n_int < 0L) {
    gdpar_abort(
      sprintf(
        paste(
          "Inconsistent B-spline configuration: df = %d must be at least",
          "degree = %d."
        ),
        W$df, W$degree
      ),
      class = "gdpar_input_error"
    )
  }
  if (n_int == 0L) {
    W$knots <- numeric(0L)
    return(W)
  }
  grid <- seq(W$boundary_knots[1L], W$boundary_knots[2L],
              length.out = n_int + 2L)
  W$knots <- grid[-c(1L, n_int + 2L)]
  W
}

#' Internal: assemble the augmented B-spline knot vector for Stan
#'
#' Sub-phase 8.3.8. Stan-side Cox--de Boor recursion consumes the
#' augmented knot vector \eqn{t} of length \eqn{n_{int} + 2(D + 1)},
#' built as \code{c(rep(b_lo, D+1), interior, rep(b_hi, D+1))} where
#' \eqn{D} is the degree, \eqn{b_{lo}/b_{hi}} the boundary knots and
#' \eqn{interior} the resolved interior knots. The total number of
#' Cox--de Boor B-splines is then \eqn{n_{int} + D + 1}, of which the
#' first one is dropped to mirror \code{splines::bs(intercept = FALSE)};
#' the returned \code{df} per coordinate equals \eqn{n_{int} + D}.
#'
#' @keywords internal
#' @noRd
.gdpar_bspline_knots_full <- function(W) {
  if (W$type != "bspline") return(numeric(0L))
  if (is.null(W$boundary_knots)) {
    gdpar_abort(
      "Cannot build augmented knot vector: 'boundary_knots' is NULL.",
      class = "gdpar_internal_error"
    )
  }
  W_resolved <- .gdpar_resolve_bspline_knots(W)
  knots_int <- if (is.null(W_resolved$knots)) numeric(0L) else W_resolved$knots
  c(rep(W_resolved$boundary_knots[1L], W_resolved$degree + 1L),
    knots_int,
    rep(W_resolved$boundary_knots[2L], W_resolved$degree + 1L))
}

#' Internal: estimate a conservative range for theta_ref from the outcome
#'
#' Sub-phase 8.3.8. The boundary-range validator (D4 of the scoping)
#' needs an a-priori range for theta_ref to compare against
#' \code{boundary_knots}. We use the canonical-link inverse of the
#' outcome as the proxy for theta_ref, expanded by one empirical
#' standard deviation on each side: identity for stan_id 1, logit for
#' stan_ids whose slot 1 has link logit (Beta, Bernoulli), log for
#' stan_ids whose slot 1 has link log (Poisson, NB, Gamma, Tweedie,
#' Lognormal, mixtures). Unknown / exotic stan_ids fall back to the
#' identity proxy with a warning.
#'
#' @keywords internal
#' @noRd
.gdpar_estimate_theta_ref_proxy_range <- function(y_vec, family) {
  z_raw <- as.numeric(y_vec)
  z_raw <- z_raw[is.finite(z_raw)]
  if (length(z_raw) < 2L) {
    return(c(-1, 1))
  }
  link_id <- tryCatch(
    .gdpar_canonical_inv_link_id_slot1(family$stan_id),
    error = function(e) 0L
  )
  z <- if (link_id == 0L) {
    z_raw
  } else if (link_id == 1L) {
    z_clamp <- pmin(pmax(z_raw, 1e-6), 1 - 1e-6)
    log(z_clamp / (1 - z_clamp))
  } else {
    log(pmax(z_raw, 1e-6))
  }
  m_rng <- range(z, na.rm = TRUE)
  s_z <- stats::sd(z, na.rm = TRUE)
  if (!is.finite(s_z) || s_z == 0) s_z <- abs(diff(m_rng)) / 4
  if (!is.finite(s_z) || s_z == 0) s_z <- 1
  c(m_rng[1L] - s_z, m_rng[2L] + s_z)
}

#' Internal: assemble the W-related Stan data fields per W basis type
#'
#' Sub-phase 8.3.8 (2026-05-22). Unified entry point used by the three
#' \code{assemble_stan_data*} variants (K = 1 + p = 1; K = 1 + p >= 1;
#' K > 1 + p = 1). Encapsulates the dispatch over \code{W$type} into a
#' single integer identifier \code{W_type_id} (0 = W off, 1 = polynomial,
#' 2 = B-spline) plus the auxiliary fields the Stan-side Cox--de Boor
#' recursion consumes (\code{W_n_knots_full}, \code{W_knots_full},
#' \code{W_degree}). For B-splines the helper also (i) resolves interior
#' knots from \code{df + boundary_knots} when knots are not explicit and
#' (ii) calls the boundary-range validator \code{.gdpar_validate_bspline_
#' boundary_range} with the projected \code{theta_ref} range derived
#' from the outcome via \code{.gdpar_estimate_theta_ref_proxy_range},
#' enforcing the contract that Stan-side relies on (D4 = delta+alpha of
#' the scoping). Returns a named list of four fields suitable for
#' \code{c()}-merging into the stan_data list.
#'
#' @keywords internal
#' @noRd
.gdpar_resolve_W_stan_data <- function(W, use_W, y_vec, family) {
  empty <- list(
    W_type_id      = 0L,
    W_n_knots_full = 0L,
    W_knots_full   = numeric(0L),
    W_degree       = 0L
  )
  if (use_W == 0L || is.null(W)) {
    return(empty)
  }
  if (identical(W$type, "polynomial")) {
    return(list(
      W_type_id      = 1L,
      W_n_knots_full = 0L,
      W_knots_full   = numeric(0L),
      W_degree       = as.integer(W$degree)
    ))
  }
  if (identical(W$type, "bspline")) {
    if (is.null(W$boundary_knots)) {
      gdpar_abort(
        paste(
          "B-spline W_basis requires 'boundary_knots' for the Stan-side",
          "fit path. Construct W with",
          "W_basis(type = 'bspline', degree = ..., (knots or df),",
          "boundary_knots = c(lo, hi))."
        ),
        class = "gdpar_input_error"
      )
    }
    W_resolved <- .gdpar_resolve_bspline_knots(W)
    proj_rng <- .gdpar_estimate_theta_ref_proxy_range(y_vec, family)
    .gdpar_validate_bspline_boundary_range(W_resolved, proj_rng)
    knots_full <- .gdpar_bspline_knots_full(W_resolved)
    return(list(
      W_type_id      = 2L,
      W_n_knots_full = length(knots_full),
      W_knots_full   = as.double(knots_full),
      W_degree       = as.integer(W_resolved$degree)
    ))
  }
  gdpar_abort(
    sprintf(
      paste(
        "Unsupported W basis type '%s' in stan_data resolver. Supported:",
        "'polynomial', 'bspline'. User-defined bases (type = 'user') are",
        "not wired to the Stan-side fit path."
      ),
      W$type
    ),
    class = "gdpar_unsupported_feature_error",
    data = list(W_type = W$type)
  )
}

#' Internal: validate that boundary_knots contains the projected theta_ref range
#'
#' Sub-phase 8.3.8, decision D4 of the scoping (delta + alpha): the
#' R-side enforces the contract "projected theta_ref range strictly
#' inside boundary_knots"; the Stan-side Cox--de Boor recursion then
#' assumes the contract and skips runtime checks. On violation, abort
#' with \code{gdpar_input_error} listing both the boundary and the
#' projected range so the user can either widen the boundary or trim
#' the data.
#'
#' @keywords internal
#' @noRd
.gdpar_validate_bspline_boundary_range <- function(W,
                                                    projected_range) {
  if (W$type != "bspline") return(invisible(NULL))
  if (is.null(W$boundary_knots)) return(invisible(NULL))
  rng <- as.numeric(projected_range)
  if (length(rng) != 2L || any(!is.finite(rng))) {
    gdpar_abort(
      "projected_range must be a length-2 finite numeric vector.",
      class = "gdpar_internal_error",
      data = list(received = projected_range)
    )
  }
  b <- W$boundary_knots
  if (rng[1L] < b[1L] || rng[2L] > b[2L]) {
    gdpar_abort(
      sprintf(
        paste(
          "B-spline boundary_knots = [%g, %g] do not strictly contain",
          "the projected theta_ref range [%g, %g]. Widen boundary_knots",
          "or restrict the outcome range to fit inside the basis support."
        ),
        b[1L], b[2L], rng[1L], rng[2L]
      ),
      class = "gdpar_input_error",
      data = list(boundary_knots = b,
                  projected_range = rng)
    )
  }
  invisible(NULL)
}

#' Start a chainable builder for an AMM specification
#'
#' Initialise an empty \code{amm_builder} object that can be incrementally
#' configured through a sequence of setter calls and finalised into an
#' \code{\link{amm_spec}} via \code{\link{as_amm_spec}}. The builder
#' pattern is intended for low-verbosity programmatic construction of
#' specifications in user scripts and serialisation routines (see
#' upcoming \code{amm_save_spec} and \code{amm_load_spec}). Functionally
#' equivalent to a direct call to \code{\link{amm_spec}}; the builder
#' adds no new modelling capability.
#'
#' @param p Positive integer giving the dimension of the per-individual
#'   parameter vector \eqn{\theta_i}. Defaults to \code{1L} (scalar path,
#'   backward compatible with the direct \code{\link{amm_spec}} call
#'   without \code{p}).
#'
#' @return An object of class \code{amm_builder} with components
#'   \code{p}, \code{dims} (a \code{\link[=dimwise]{dims_spec}} object
#'   initialised with NULL base and no overrides), \code{W} (NULL), and
#'   \code{x_vars} (NULL). The object is intended to be passed through a
#'   chain of \code{amm_set_*} calls and terminated by
#'   \code{\link{as_amm_spec}}.
#'
#' @details
#' The builder is a thin wrapper around the existing \code{dims_spec}
#' semantics: every \code{amm_set_a_uniform} / \code{amm_set_b_uniform}
#' call mutates the base of the embedded \code{dims_spec}, and every
#' \code{amm_set_a} / \code{amm_set_b} call layers a per-dimension
#' override on top. Overrides survive subsequent uniform changes, in
#' direct correspondence to \code{\link{dimwise}} composed with
#' \code{\link{override}}.
#'
#' Construction is bifurcated by \code{p} at the finalisation step, not
#' at builder time: when \code{as_amm_spec} is called with \code{p = 1L},
#' the embedded \code{dims_spec} is resolved to a length-one list and
#' the resulting \code{a} and \code{b} formulas are passed to
#' \code{\link{amm_spec}} on the scalar path; with \code{p > 1L} the
#' embedded \code{dims_spec} is passed directly to the multivariate
#' path. This keeps a single source of truth for per-dimension
#' semantics.
#'
#' @section Methodological notes:
#' The builder exposes the same structural asymmetry between
#' per-dimension components (\code{a}, \code{b}) and the cross-dimension
#' component (\code{W}) as \code{\link{amm_spec}}: dedicated
#' per-dimension setters are provided for \code{a} and \code{b}, while
#' \code{W} is set globally via \code{\link{amm_set_W}}. Declaring
#' \code{W} per-dimension would silently restrict the model class to
#' the separable sub-class, which is rejected by construction.
#'
#' @examples
#' spec <- amm_build(p = 1L) |>
#'   amm_set_a_uniform(~ x1 + x2) |>
#'   amm_set_b_uniform(~ x1) |>
#'   amm_set_W(W_basis(type = "polynomial", degree = 2)) |>
#'   as_amm_spec()
#' print(spec)
#'
#' spec_mv <- amm_build(p = 3L) |>
#'   amm_set_a_uniform(~ x1 + x2) |>
#'   amm_set_a(k = 2L, ~ x1) |>
#'   amm_set_b_uniform(~ x1) |>
#'   as_amm_spec()
#' print(spec_mv)
#'
#' @seealso \code{\link{amm_set_a_uniform}}, \code{\link{amm_set_b_uniform}},
#'   \code{\link{amm_set_a}}, \code{\link{amm_set_b}},
#'   \code{\link{amm_set_W}}, \code{\link{amm_set_x_vars}},
#'   \code{\link{as_amm_spec}}, \code{\link{amm_spec}},
#'   \code{\link{dimwise}}, \code{\link{override}}
#' @export
amm_build <- function(p = 1L) {
  assert_count(p, "p")
  obj <- list(
    p = as.integer(p),
    dims = dimwise(a = NULL, b = NULL),
    W = NULL,
    x_vars = NULL
  )
  class(obj) <- c("amm_builder", "list")
  obj
}

#' Set the uniform additive basis on an AMM builder
#'
#' Replace the base additive formula of the embedded \code{dims_spec} with
#' the supplied one-sided formula (or \code{NULL} to disable the additive
#' component on the base). Per-dimension overrides previously registered
#' via \code{\link{amm_set_a}} are preserved and continue to take
#' precedence at their respective indices.
#'
#' @param builder An object of class \code{amm_builder}.
#' @param a One-sided formula or \code{NULL}. Becomes the new base of the
#'   additive component, applied uniformly to every dimension that lacks
#'   an explicit override.
#'
#' @return The modified \code{amm_builder}, returned invisibly for the
#'   pipe convention but suitable for direct inspection.
#'
#' @examples
#' b <- amm_build(p = 2L) |>
#'   amm_set_a(k = 2L, ~ x1) |>
#'   amm_set_a_uniform(~ x1 + x2)
#' print(b)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{amm_set_a}}
#' @export
amm_set_a_uniform <- function(builder, a) {
  assert_inherits(builder, "amm_builder", "builder")
  assert_one_sided_formula(a, "a", allow_null = TRUE)
  builder$dims$base$a <- a
  builder
}

#' Set the uniform multiplicative basis on an AMM builder
#'
#' Replace the base multiplicative formula of the embedded
#' \code{dims_spec} with the supplied one-sided formula (or \code{NULL}
#' to disable the multiplicative component on the base). Per-dimension
#' overrides previously registered via \code{\link{amm_set_b}} are
#' preserved.
#'
#' @param builder An object of class \code{amm_builder}.
#' @param b One-sided formula or \code{NULL}. Becomes the new base of the
#'   multiplicative component, applied uniformly to every dimension that
#'   lacks an explicit override.
#'
#' @return The modified \code{amm_builder}.
#'
#' @examples
#' b <- amm_build(p = 2L) |>
#'   amm_set_b_uniform(~ x1)
#' print(b)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{amm_set_b}}
#' @export
amm_set_b_uniform <- function(builder, b) {
  assert_inherits(builder, "amm_builder", "builder")
  assert_one_sided_formula(b, "b", allow_null = TRUE)
  builder$dims$base$b <- b
  builder
}

#' Set a per-dimension additive basis override on an AMM builder
#'
#' Register a per-dimension override of the additive component for index
#' \code{k}. The override replaces the uniform base (set via
#' \code{\link{amm_set_a_uniform}}, or NULL if never set) for dimension
#' \code{k} only. Calling \code{amm_set_a} twice with the same \code{k}
#' replaces the previous override.
#'
#' @param builder An object of class \code{amm_builder}.
#' @param k Positive integer in \code{1:p}. Indexes the dimension to
#'   override.
#' @param a One-sided formula or \code{NULL}. The latter disables the
#'   additive component for dimension \code{k} only.
#'
#' @return The modified \code{amm_builder}.
#'
#' @examples
#' b <- amm_build(p = 3L) |>
#'   amm_set_a_uniform(~ x1 + x2) |>
#'   amm_set_a(k = 2L, ~ x1) |>
#'   amm_set_a(k = 3L, NULL)
#' print(b)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{amm_set_a_uniform}},
#'   \code{\link{override}}
#' @export
amm_set_a <- function(builder, k, a) {
  assert_inherits(builder, "amm_builder", "builder")
  assert_count(k, "k", max = builder$p)
  assert_one_sided_formula(a, "a", allow_null = TRUE)
  builder$dims <- override(builder$dims, k = k, a = a)
  builder
}

#' Set a per-dimension multiplicative basis override on an AMM builder
#'
#' Register a per-dimension override of the multiplicative component for
#' index \code{k}. The override replaces the uniform base (set via
#' \code{\link{amm_set_b_uniform}}, or NULL if never set) for dimension
#' \code{k} only. Calling \code{amm_set_b} twice with the same \code{k}
#' replaces the previous override.
#'
#' @param builder An object of class \code{amm_builder}.
#' @param k Positive integer in \code{1:p}.
#' @param b One-sided formula or \code{NULL}.
#'
#' @return The modified \code{amm_builder}.
#'
#' @examples
#' b <- amm_build(p = 2L) |>
#'   amm_set_b_uniform(~ x1) |>
#'   amm_set_b(k = 2L, NULL)
#' print(b)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{amm_set_b_uniform}},
#'   \code{\link{override}}
#' @export
amm_set_b <- function(builder, k, b) {
  assert_inherits(builder, "amm_builder", "builder")
  assert_count(k, "k", max = builder$p)
  assert_one_sided_formula(b, "b", allow_null = TRUE)
  builder$dims <- override(builder$dims, k = k, b = b)
  builder
}

#' Set the modulating basis on an AMM builder
#'
#' Store a \code{\link{W_basis}} object as the modulating component of
#' the specification under construction, or clear it by passing
#' \code{NULL}. The modulating component is global to all dimensions of
#' \eqn{\theta_i} and is therefore stored as a single top-level slot of
#' the builder.
#'
#' @param builder An object of class \code{amm_builder}.
#' @param W A \code{W_basis} object, or \code{NULL} to disable the
#'   modulating component.
#'
#' @return The modified \code{amm_builder}.
#'
#' @examples
#' wb <- W_basis(type = "polynomial", degree = 2)
#' b <- amm_build(p = 1L) |>
#'   amm_set_a_uniform(~ x1) |>
#'   amm_set_W(wb)
#' print(b)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{W_basis}}
#' @export
amm_set_W <- function(builder, W) {
  assert_inherits(builder, "amm_builder", "builder")
  if (!is.null(W)) {
    assert_inherits(W, "W_basis", "W")
  }
  builder$W <- W
  builder
}

#' Set the covariate names used by the modulating component on an AMM builder
#'
#' Record the character vector identifying the covariates that enter the
#' modulating component as the linear factor \code{x} in
#' \eqn{W(\theta) x}, or clear it by passing \code{NULL}. The value is
#' forwarded to \code{\link{amm_spec}} at finalisation time. When NULL,
#' the package uses the covariates derived from the right-hand side of
#' the model formula passed to \code{\link{gdpar}}.
#'
#' @param builder An object of class \code{amm_builder}.
#' @param x_vars Character vector with the names of the covariates, or
#'   \code{NULL}.
#'
#' @return The modified \code{amm_builder}.
#'
#' @examples
#' b <- amm_build(p = 1L) |>
#'   amm_set_a_uniform(~ x1 + x2) |>
#'   amm_set_x_vars(c("x1", "x2"))
#' print(b)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{amm_spec}}
#' @export
amm_set_x_vars <- function(builder, x_vars) {
  assert_inherits(builder, "amm_builder", "builder")
  if (!is.null(x_vars)) {
    if (!is.character(x_vars) || length(x_vars) < 1L) {
      gdpar_abort(
        "Argument 'x_vars' must be a non-empty character vector or NULL.",
        class = "gdpar_input_error",
        data = list(argument = "x_vars", received = x_vars)
      )
    }
  }
  builder$x_vars <- x_vars
  builder
}

#' Finalise an AMM builder into an amm_spec
#'
#' Validate the accumulated state of an \code{amm_builder} and convert it
#' to an \code{\link{amm_spec}} object. Dispatches on \code{p}: the
#' scalar path (\code{p = 1L}) resolves the embedded \code{dims_spec} to
#' the single per-dimension entry and forwards \code{a}, \code{b},
#' \code{W}, \code{x_vars} to \code{\link{amm_spec}}; the multivariate
#' path (\code{p > 1L}) forwards the embedded \code{dims_spec} directly.
#' All structural validation (range of override indices, consistency
#' between \code{p} and the supplied components) is delegated to
#' \code{\link{amm_spec}}.
#'
#' @param builder An object of class \code{amm_builder}.
#'
#' @return An object of class \code{amm_spec}; see \code{\link{amm_spec}}
#'   for the slot layout.
#'
#' @examples
#' spec <- amm_build(p = 2L) |>
#'   amm_set_a_uniform(~ x1 + x2) |>
#'   amm_set_b_uniform(~ x1) |>
#'   as_amm_spec()
#' print(spec)
#'
#' @seealso \code{\link{amm_build}}, \code{\link{amm_spec}}
#' @export
as_amm_spec <- function(builder) {
  assert_inherits(builder, "amm_builder", "builder")
  p <- builder$p
  if (p == 1L) {
    resolved <- resolve_dims_spec(builder$dims, 1L)
    amm_spec(
      a = resolved[[1L]][["a"]],
      b = resolved[[1L]][["b"]],
      W = builder$W,
      x_vars = builder$x_vars,
      p = 1L
    )
  } else {
    amm_spec(
      W = builder$W,
      x_vars = builder$x_vars,
      p = p,
      dims = builder$dims
    )
  }
}

#' Print method for amm_builder objects
#'
#' @param x An object of class \code{amm_builder}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.amm_builder <- function(x, ...) {
  cat("<amm_builder> p = ", x$p, "\n", sep = "")
  cat("  a (base)    : ",
      if (is.null(x$dims$base$a)) "NULL" else deparse(x$dims$base$a),
      "\n", sep = "")
  cat("  b (base)    : ",
      if (is.null(x$dims$base$b)) "NULL" else deparse(x$dims$base$b),
      "\n", sep = "")
  if (length(x$dims$overrides) > 0L) {
    cat("  overrides   :\n")
    keys <- sort(as.integer(names(x$dims$overrides)))
    for (k in keys) {
      ov <- x$dims$overrides[[as.character(k)]]
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
    cat("  overrides   : <none>\n")
  }
  cat("  W           : ",
      if (is.null(x$W)) "NULL" else paste0("W_basis(type = '",
                                            x$W$type, "')"),
      "\n", sep = "")
  if (!is.null(x$x_vars)) {
    cat("  x_vars      : ",
        paste(x$x_vars, collapse = ", "), "\n", sep = "")
  } else {
    cat("  x_vars      : NULL (inherited from gdpar() formula)\n")
  }
  invisible(x)
}

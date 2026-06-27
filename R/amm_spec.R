#' Specify the AMM canonical decomposition
#'
#' Declare which components of the Additive-Multiplicative-Modulated
#' canonical form
#' \eqn{\Delta(x, \theta) = a(x) + b(x) \odot \theta + W(\theta) x}
#' are active and how they are parametrized. The result is consumed by
#' \code{\link{gdpar}} to assemble the design matrices and the Stan
#' model.
#'
#' Two construction paths exist, selected by the dimension \code{p} of
#' the per-individual parameter vector \eqn{\theta_i}:
#' the scalar path (\code{p = 1L}, the default) uses \code{a} and \code{b}
#' as one-sided formulas directly; the multivariate path (\code{p > 1L})
#' uses the \code{dims} argument together with the \code{\link{dimwise}}
#' and \code{\link{override}} helpers to declare per-dimension
#' specifications. The two paths are deliberately bifurcated to keep
#' contract semantics clean and backward compatibility absolute.
#'
#' @param a Scalar path only (\code{p = 1L}). One-sided formula declaring
#'   the basis of the additive component on the covariate space,
#'   evaluated via \code{\link[stats]{model.matrix}}. Use \code{NULL} to
#'   disable the additive component (degenerate Level 0 or Level 1 with
#'   only the modulating term). Must be \code{NULL} when \code{p > 1L}.
#' @param b Scalar path only (\code{p = 1L}). One-sided formula declaring
#'   the basis of the multiplicative component on the covariate space.
#'   Must be \code{NULL} when \code{p > 1L}.
#' @param W An object of class \code{W_basis} declaring the basis of
#'   the modulating component on the parameter space, or \code{NULL} to
#'   disable the modulating component. See \code{\link{W_basis}}. The
#'   modulating component couples all dimensions of
#'   \eqn{\theta_{\mathrm{ref}}} in the canonical form and is therefore
#'   declared as a single basis object regardless of \code{p}.
#' @param x_vars Character vector with the names of the covariates that
#'   enter the modulating component as the linear factor x in
#'   \eqn{W(\theta) x}. When \code{NULL} (default), the package uses the
#'   same covariates as those entering the right-hand side of the model
#'   formula passed to \code{\link{gdpar}}.
#' @param p Positive integer giving the dimension of the per-individual
#'   parameter vector \eqn{\theta_i}. Defaults to \code{1L} (scalar path,
#'   backward compatible with all previous specifications).
#' @param dims Multivariate path only (\code{p > 1L}). Either a
#'   \code{dims_spec} object produced by \code{\link{dimwise}} (possibly
#'   composed with \code{\link{override}}), or a plain list of length
#'   \code{p} where each entry is itself a list with components \code{a}
#'   and \code{b} (each a one-sided formula or \code{NULL}). Must be
#'   \code{NULL} when \code{p == 1L}. Bare formulas are rejected, to
#'   prevent silent recycling.
#'
#' @return An object of class \code{amm_spec} with components \code{a},
#'   \code{b}, \code{W}, \code{x_vars}, \code{level}, \code{p}, and
#'   \code{dims}. On the scalar path (\code{p = 1L}), \code{dims} is
#'   \code{NULL} and \code{a}, \code{b} hold the formulas. On the
#'   multivariate path (\code{p > 1L}), \code{a} and \code{b} are
#'   \code{NULL} and \code{dims} holds the resolved per-dimension list of
#'   length \code{p}, each entry being a list with \code{a} and \code{b}.
#'
#' @details
#' The additive and multiplicative bases are declared as one-sided R
#' formulas; the package uses \code{\link[stats]{model.matrix}} to build
#' the design matrices from the data frame supplied to \code{\link{gdpar}}.
#' Each formula is evaluated without an intercept column; the centering
#' conditions (C2) and (C3) of Block 1 are enforced empirically by
#' column-wise centering of the design matrices Z_a and Z_b inside the
#' constructor (each column has its sample mean subtracted, so that
#' colMeans(Z_a) = 0 and colMeans(Z_b) = 0 exactly). Because
#' E_mu[a(X)] = colMeans(Z_a) * a_coef under the empirical mu used by
#' Path 1, the centering of Z_a alone satisfies (C2) for any choice of
#' a_coef in the full J_a-dimensional Euclidean space; analogously for
#' (C3). No additional restriction on the basis coefficients is imposed.
#'
#' On the multivariate path, the additive and multiplicative components
#' factor per dimension of \eqn{\theta_i} because they depend only on the
#' covariates \eqn{x_i} (and therefore on coordinates of
#' \eqn{\mathcal{X}}, not of \eqn{\Theta}). The modulating component
#' \eqn{W} depends on the full vector \eqn{\theta_{\mathrm{ref}}} and
#' therefore couples the dimensions; declaring \eqn{W} per-dimension
#' would silently restrict the model class to the separable sub-class.
#' The package therefore enforces the structural asymmetry between
#' (\code{a}, \code{b}) and \code{W} at the API level: \code{dims} for
#' the former, a single top-level \code{W} argument for the latter.
#'
#' @section Methodological notes:
#' The AMM level implied by the specification is recorded on the
#' returned object: Level 0 if every component is NULL across all
#' dimensions, Level 1 if only the additive component is active in any
#' dimension and no multiplicative or modulating term is active anywhere,
#' and Level 2 otherwise. The identifiability theorems of Block 1 apply
#' at each level under the linearity assumption (LIN), which holds
#' automatically for formula-based bases (linear subspaces of L^2_0(\\mu))
#' and for the polynomial and B-spline cases of \code{\link{W_basis}}.
#'
#' Cross-component non-identifiability is not detected by this
#' constructor; it is detected at the Gram-matrix diagnostic of
#' \code{\link{gdpar_check_identifiability}}, which is called
#' automatically before fitting. For \code{p > 1L}, an additional
#' cross-dimension identifiability condition applies (see the planned
#' vignette \code{vop02_arbitrary_p}).
#'
#' @section Dependencies:
#' Uses \code{\link[stats]{terms}} for parsing the formulas and
#' \code{\link[stats]{model.matrix}} downstream for building the design
#' matrices.
#'
#' @references
#' See \code{vignette("v01_amm_identifiability", package = "gdpar")},
#' Section 3 (the AMM hierarchy) and Section 5 (standing assumptions
#' (C1)-(C6)).
#'
#' @examples
#' spec <- amm_spec(
#'   a = ~ x1 + x2,
#'   b = ~ x1,
#'   W = W_basis(type = "polynomial", degree = 2)
#' )
#' print(spec)
#'
#' spec_mv <- amm_spec(
#'   p    = 2L,
#'   dims = dimwise(a = ~ x1 + x2, b = ~ x1),
#'   W    = W_basis(type = "polynomial", degree = 2)
#' )
#' print(spec_mv)
#'
#' @seealso \code{\link{W_basis}}, \code{\link{dimwise}},
#'   \code{\link{override}}, \code{\link{gdpar_check_identifiability}},
#'   \code{\link{gdpar}}
#' @export
amm_spec <- function(a = NULL, b = NULL, W = NULL, x_vars = NULL,
                     p = 1L, dims = NULL) {
  assert_count(p, "p")
  p <- as.integer(p)
  if (!is.null(W)) {
    assert_inherits(W, "W_basis", "W")
  }
  if (!is.null(x_vars)) {
    if (!is.character(x_vars) || length(x_vars) < 1L) {
      gdpar_abort(
        "Argument 'x_vars' must be a character vector of variable names.",
        class = "gdpar_input_error"
      )
    }
  }
  if (p == 1L) {
    if (!is.null(dims)) {
      gdpar_abort(
        paste(
          "When p = 1L (default), the 'dims' argument must be NULL.",
          "Use 'a' and 'b' directly."
        ),
        class = "gdpar_input_error"
      )
    }
    assert_one_sided_formula(a, "a", allow_null = TRUE)
    assert_one_sided_formula(b, "b", allow_null = TRUE)
    level <- if (is.null(a) && is.null(b) && is.null(W)) {
      0L
    } else if (is.null(b) && is.null(W)) {
      1L
    } else {
      2L
    }
    obj <- list(
      a = a, b = b, W = W,
      x_vars = x_vars,
      level = level,
      p = 1L,
      dims = NULL
    )
    class(obj) <- c("amm_spec", "list")
    return(obj)
  }
  if (!is.null(a) || !is.null(b)) {
    gdpar_abort(
      paste(
        "When p > 1L, per-dimension specifications must be supplied via",
        "'dims' (see dimwise()). Arguments 'a' and 'b' must be NULL in",
        "the multivariate case."
      ),
      class = "gdpar_input_error"
    )
  }
  if (is.null(dims)) {
    gdpar_abort(
      paste(
        "When p > 1L, the 'dims' argument is required.",
        "See dimwise() to build it from a uniform template, optionally",
        "composed with override() for per-dimension changes."
      ),
      class = "gdpar_input_error"
    )
  }
  if (inherits(dims, "formula")) {
    gdpar_abort(
      paste(
        "A bare formula passed to 'dims' is ambiguous when p > 1L.",
        "Use dimwise() to broadcast a uniform formula to all dimensions,",
        "or supply a list of length p with per-dimension specifications."
      ),
      class = "gdpar_input_error"
    )
  }
  if (inherits(dims, "dims_spec")) {
    resolved <- resolve_dims_spec(dims, p)
  } else if (is.list(dims)) {
    if (length(dims) != p) {
      gdpar_abort(
        sprintf(
          paste(
            "When 'dims' is a plain list and p = %d, the list must have",
            "length %d; received length %d."
          ),
          p, p, length(dims)
        ),
        class = "gdpar_input_error",
        data = list(p = p, received_length = length(dims))
      )
    }
    resolved <- vector("list", p)
    for (k in seq_len(p)) {
      entry <- dims[[k]]
      if (!is.list(entry)) {
        gdpar_abort(
          sprintf(
            paste(
              "Entry %d of 'dims' must be a list with components 'a' and",
              "'b'; received an object of class %s."
            ),
            k, paste(sQuote(class(entry)), collapse = ", ")
          ),
          class = "gdpar_input_error"
        )
      }
      a_k <- if ("a" %in% names(entry)) entry[["a"]] else NULL
      b_k <- if ("b" %in% names(entry)) entry[["b"]] else NULL
      assert_one_sided_formula(
        a_k, sprintf("dims[[%d]]$a", k), allow_null = TRUE
      )
      assert_one_sided_formula(
        b_k, sprintf("dims[[%d]]$b", k), allow_null = TRUE
      )
      resolved[[k]] <- list(a = a_k, b = b_k)
    }
  } else {
    gdpar_abort(
      sprintf(
        paste(
          "Argument 'dims' must be a dims_spec (from dimwise()) or a list",
          "of length %d; received object of class %s."
        ),
        p, paste(sQuote(class(dims)), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(p = p, received_class = class(dims))
    )
  }
  any_a <- any(vapply(resolved, function(d) !is.null(d$a), logical(1L)))
  any_b <- any(vapply(resolved, function(d) !is.null(d$b), logical(1L)))
  level <- if (!any_a && !any_b && is.null(W)) {
    0L
  } else if (!any_b && is.null(W)) {
    1L
  } else {
    2L
  }
  obj <- list(
    a = NULL, b = NULL, W = W,
    x_vars = x_vars,
    level = level,
    p = p,
    dims = resolved
  )
  class(obj) <- c("amm_spec", "list")
  obj
}

#' Internal: build the design matrices for an AMM specification
#'
#' Materializes the design matrices implied by an \code{amm_spec}
#' object and a data frame: the additive basis Z_a, the multiplicative
#' basis Z_b, and the linear factor X used by the modulating
#' component. Covariates are centered prior to matrix construction
#' to enforce assumption (C1) of Block 1 at the empirical level.
#'
#' @param amm An object of class \code{amm_spec}.
#' @param data A data frame containing the variables referenced by
#'   \code{amm}.
#' @param formula_rhs A formula or character vector identifying the
#'   covariate columns of \code{data} used as the linear factor x
#'   when \code{amm$x_vars} is NULL.
#'
#' @return A list with components \code{Z_a} (numeric matrix, n by
#'   J_a), \code{Z_b} (numeric matrix, n by J_b), \code{X} (numeric
#'   matrix, n by d), \code{Z_a_means}, \code{Z_b_means},
#'   \code{X_means}, \code{X_sds}, plus the column names of each
#'   design block. Centering uses the column means; the X block is
#'   additionally scaled to unit standard deviation per column. The
#'   centering is recorded so that predictions on new data can apply
#'   the same transformation.
#' @keywords internal
#' @noRd
build_amm_design <- function(amm, data, formula_rhs) {
  assert_inherits(amm, "amm_spec", "amm")
  assert_data_frame(data, "data")
  if (!is.null(amm$p) && amm$p > 1L) {
    return(.build_amm_design_multi(amm, data, formula_rhs))
  }
  needed_vars <- character(0)
  if (!is.null(amm$a)) {
    needed_vars <- c(needed_vars, all.vars(amm$a))
  }
  if (!is.null(amm$b)) {
    needed_vars <- c(needed_vars, all.vars(amm$b))
  }
  if (!is.null(amm$x_vars)) {
    needed_vars <- c(needed_vars, amm$x_vars)
  }
  needed_vars <- unique(needed_vars)
  needed_vars <- intersect(needed_vars, colnames(data))
  if (length(needed_vars) > 0L) {
    has_na <- vapply(needed_vars,
                     function(v) any(is.na(data[[v]])),
                     logical(1L))
    if (any(has_na)) {
      gdpar_abort(
        sprintf(
          "Covariate(s) used by the AMM specification contain missing values: %s. Path 1 does not impute.",
          paste(sQuote(needed_vars[has_na]), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
  }
  n <- nrow(data)
  if (!is.null(amm$a)) {
    Z_a_full <- stats::model.matrix(stats::update(amm$a, ~ . + 0), data)
    Z_a_full <- as.matrix(Z_a_full)
    Z_a_means <- colMeans(Z_a_full)
    Z_a <- sweep(Z_a_full, 2L, Z_a_means, "-")
  } else {
    Z_a <- matrix(numeric(0), nrow = n, ncol = 0)
    Z_a_means <- numeric(0)
  }
  if (!is.null(amm$b)) {
    Z_b_full <- stats::model.matrix(stats::update(amm$b, ~ . + 0), data)
    Z_b_full <- as.matrix(Z_b_full)
    Z_b_means <- colMeans(Z_b_full)
    Z_b <- sweep(Z_b_full, 2L, Z_b_means, "-")
  } else {
    Z_b <- matrix(numeric(0), nrow = n, ncol = 0)
    Z_b_means <- numeric(0)
  }
  if (!is.null(amm$W)) {
    x_vars <- if (!is.null(amm$x_vars)) {
      amm$x_vars
    } else if (inherits(formula_rhs, "formula")) {
      attr(stats::terms(formula_rhs, data = data), "term.labels")
    } else if (is.character(formula_rhs)) {
      formula_rhs
    } else {
      character(0)
    }
    if (length(x_vars) == 0L) {
      gdpar_abort(
        paste(
          "The modulating component W is active but no covariates were",
          "identified for the linear factor x. Either supply 'x_vars'",
          "to amm_spec or include covariates in the model formula."
        ),
        class = "gdpar_input_error"
      )
    }
    missing_x <- setdiff(x_vars, colnames(data))
    if (length(missing_x) > 0L) {
      gdpar_abort(
        sprintf(
          "Variables required by the modulating component are missing in data: %s.",
          paste(sQuote(missing_x), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
    X_full <- as.matrix(data[, x_vars, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_means <- colMeans(X_full)
    X_centered <- sweep(X_full, 2L, X_means, "-")
    X_sds <- apply(X_centered, 2L, stats::sd)
    if (any(X_sds == 0)) {
      bad <- x_vars[X_sds == 0]
      gdpar_abort(
        sprintf(
          "Constant covariate(s) cannot enter the modulating component: %s.",
          paste(sQuote(bad), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
    X <- sweep(X_centered, 2L, X_sds, "/")
  } else {
    X <- matrix(numeric(0), nrow = n, ncol = 0)
    X_means <- numeric(0)
    X_sds <- numeric(0)
    x_vars <- character(0)
  }
  list(
    Z_a = Z_a, Z_b = Z_b, X = X,
    Z_a_means = Z_a_means, Z_b_means = Z_b_means,
    X_means = X_means, X_sds = X_sds,
    Z_a_names = colnames(Z_a), Z_b_names = colnames(Z_b),
    X_names = if (length(x_vars) > 0L) x_vars else character(0)
  )
}

#' Internal: build the design matrices for an AMM specification with p > 1
#'
#' Multivariate path of \code{build_amm_design}. Iterates over the
#' per-coordinate specifications in \code{amm$dims} and builds one
#' centered design matrix \code{Z_a_k} and \code{Z_b_k} per coordinate
#' \eqn{k = 1, \dots, p}. The modulating component \code{X} is shared
#' across coordinates because \eqn{W(\theta_{\mathrm{ref}}) x} depends
#' only on the global covariate vector \eqn{x}, not on the coordinate
#' index.
#'
#' The per-coordinate design matrices are ragged in their column counts
#' (\eqn{J_{a,k}} may differ across \eqn{k}). The padded layout
#' required by the Stan template is produced by
#' \code{assemble_stan_data} from these ragged lists.
#'
#' @param amm An object of class \code{amm_spec} with \code{p > 1L}.
#' @param data A data frame containing the variables referenced by
#'   \code{amm}.
#' @param formula_rhs A formula or character vector identifying the
#'   covariate columns of \code{data} used as the linear factor x
#'   when \code{amm$x_vars} is NULL.
#' @return A list with components \code{p}, \code{Z_a_list},
#'   \code{Z_b_list} (each of length \code{p}), \code{X}, and the
#'   per-coordinate means and names lists, plus shared \code{X_means},
#'   \code{X_sds}, \code{X_names}.
#' @keywords internal
#' @noRd
.build_amm_design_multi <- function(amm, data, formula_rhs) {
  p <- amm$p
  n <- nrow(data)
  needed_vars <- character(0)
  for (k in seq_len(p)) {
    a_k <- amm$dims[[k]]$a
    b_k <- amm$dims[[k]]$b
    if (!is.null(a_k)) {
      needed_vars <- c(needed_vars, all.vars(a_k))
    }
    if (!is.null(b_k)) {
      needed_vars <- c(needed_vars, all.vars(b_k))
    }
  }
  if (!is.null(amm$x_vars)) {
    needed_vars <- c(needed_vars, amm$x_vars)
  }
  needed_vars <- unique(needed_vars)
  needed_vars <- intersect(needed_vars, colnames(data))
  if (length(needed_vars) > 0L) {
    has_na <- vapply(needed_vars,
                     function(v) any(is.na(data[[v]])),
                     logical(1L))
    if (any(has_na)) {
      gdpar_abort(
        sprintf(
          "Covariate(s) used by the AMM specification contain missing values: %s. Path 1 does not impute.",
          paste(sQuote(needed_vars[has_na]), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
  }
  Z_a_list <- vector("list", p)
  Z_b_list <- vector("list", p)
  Z_a_means_list <- vector("list", p)
  Z_b_means_list <- vector("list", p)
  Z_a_names_list <- vector("list", p)
  Z_b_names_list <- vector("list", p)
  for (k in seq_len(p)) {
    a_k <- amm$dims[[k]]$a
    b_k <- amm$dims[[k]]$b
    if (!is.null(a_k)) {
      Z_full <- stats::model.matrix(stats::update(a_k, ~ . + 0), data)
      Z_full <- as.matrix(Z_full)
      means <- colMeans(Z_full)
      Z_a_list[[k]] <- sweep(Z_full, 2L, means, "-")
      Z_a_means_list[[k]] <- means
      Z_a_names_list[[k]] <- colnames(Z_full)
    } else {
      Z_a_list[[k]] <- matrix(numeric(0), nrow = n, ncol = 0L)
      Z_a_means_list[[k]] <- numeric(0)
      Z_a_names_list[[k]] <- character(0)
    }
    if (!is.null(b_k)) {
      Z_full <- stats::model.matrix(stats::update(b_k, ~ . + 0), data)
      Z_full <- as.matrix(Z_full)
      means <- colMeans(Z_full)
      Z_b_list[[k]] <- sweep(Z_full, 2L, means, "-")
      Z_b_means_list[[k]] <- means
      Z_b_names_list[[k]] <- colnames(Z_full)
    } else {
      Z_b_list[[k]] <- matrix(numeric(0), nrow = n, ncol = 0L)
      Z_b_means_list[[k]] <- numeric(0)
      Z_b_names_list[[k]] <- character(0)
    }
  }
  if (!is.null(amm$W)) {
    x_vars <- if (!is.null(amm$x_vars)) {
      amm$x_vars
    } else if (inherits(formula_rhs, "formula")) {
      attr(stats::terms(formula_rhs, data = data), "term.labels")
    } else if (is.character(formula_rhs)) {
      formula_rhs
    } else {
      character(0)
    }
    if (length(x_vars) == 0L) {
      gdpar_abort(
        paste(
          "The modulating component W is active but no covariates were",
          "identified for the linear factor x. Either supply 'x_vars'",
          "to amm_spec or include covariates in the model formula."
        ),
        class = "gdpar_input_error"
      )
    }
    missing_x <- setdiff(x_vars, colnames(data))
    if (length(missing_x) > 0L) {
      gdpar_abort(
        sprintf(
          "Variables required by the modulating component are missing in data: %s.",
          paste(sQuote(missing_x), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
    X_full <- as.matrix(data[, x_vars, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_means <- colMeans(X_full)
    X_centered <- sweep(X_full, 2L, X_means, "-")
    X_sds <- apply(X_centered, 2L, stats::sd)
    if (any(X_sds == 0)) {
      bad <- x_vars[X_sds == 0]
      gdpar_abort(
        sprintf(
          "Constant covariate(s) cannot enter the modulating component: %s.",
          paste(sQuote(bad), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
    X <- sweep(X_centered, 2L, X_sds, "/")
  } else {
    X <- matrix(numeric(0), nrow = n, ncol = 0L)
    X_means <- numeric(0)
    X_sds <- numeric(0)
    x_vars <- character(0)
  }
  list(
    p              = p,
    Z_a_list       = Z_a_list,
    Z_b_list       = Z_b_list,
    X              = X,
    Z_a_means_list = Z_a_means_list,
    Z_b_means_list = Z_b_means_list,
    X_means        = X_means,
    X_sds          = X_sds,
    Z_a_names_list = Z_a_names_list,
    Z_b_names_list = Z_b_names_list,
    X_names        = if (length(x_vars) > 0L) x_vars else character(0)
  )
}

#' Internal: build the design matrices for a K-individual specification
#'
#' Multi-parametric counterpart of \code{build_amm_design()} for the
#' regime \eqn{K > 1} with \eqn{p = 1}. Consumes the canonical named
#' list of per-slot \code{amm_spec} objects returned by
#' \code{.gdpar_formula_set_to_amm_spec_list()} (each entry corresponds
#' to one K-individual structural parameter, with \code{p = 1}) and
#' produces one centered design matrix \code{Z_a_k} and \code{Z_b_k}
#' per slot \eqn{k = 1, \dots, K}. The modulating block \code{X} is
#' shared across slots because the scope of \eqn{W} is global by
#' design decision (Sub-decision 2 of Unit 3 paso 3, 2026-05-20); the
#' x-vars span the union of covariates that appear in any slot's
#' \code{a()} or \code{b()} component, mirroring the union-of-vars
#' construction of the K = 1 path in \code{gdpar()}.
#'
#' @param amm_list_canonical Named list of \code{amm_spec} objects of
#'   length \eqn{K \ge 2}. Each entry must satisfy \code{p = 1} (or
#'   \code{NULL}, treated as \code{p = 1}); the slot name carries the
#'   canonical parameter name (e.g. \code{"mu"}, \code{"sigma"},
#'   \code{"phi"}).
#' @param data A data frame containing the variables referenced by the
#'   per-slot AMM specifications.
#' @param formula_rhs A formula or character vector identifying the
#'   covariate columns of \code{data} used as the linear factor x when
#'   the union-of-vars construction yields an empty set and W is active
#'   on at least one slot. In the canonical pathway (formula set with
#'   AMM wrappers) this argument is unused because the union of vars
#'   in the wrappers is non-empty whenever W appears.
#'
#' @return A list with components \code{K}, \code{Z_a_k_list},
#'   \code{Z_b_k_list} (each of length \code{K}), \code{X} (shared
#'   matrix), \code{Z_a_k_means_list}, \code{Z_b_k_means_list},
#'   \code{X_means}, \code{X_sds}, \code{Z_a_k_names_list},
#'   \code{Z_b_k_names_list}, \code{X_names}, plus \code{slot_names}
#'   (the names of the K slots in order).
#' @keywords internal
#' @noRd
.build_amm_design_K <- function(amm_list_canonical, data, formula_rhs) {
  if (!is.list(amm_list_canonical) || length(amm_list_canonical) < 2L) {
    gdpar_abort(
      "Internal error: .build_amm_design_K expected a named list of length K >= 2.",
      class = "gdpar_internal_error",
      data = list(received_length = length(amm_list_canonical))
    )
  }
  slot_names <- names(amm_list_canonical)
  if (is.null(slot_names) || any(!nzchar(slot_names))) {
    gdpar_abort(
      "Internal error: .build_amm_design_K expected every slot to carry a non-empty name.",
      class = "gdpar_internal_error"
    )
  }
  K <- length(amm_list_canonical)
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]
    assert_inherits(a_k, "amm_spec", paste0("amm_list_canonical[['",
                                             slot_names[k], "']]"))
    if (!is.null(a_k$p) && a_k$p > 1L) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' of the K-individual list carries amm_spec(p = %d).",
            "The K > 1 + p > 1 regime is queued for Block 8.1; the current",
            "path supports K > 1 with p = 1 per slot."
          ),
          slot_names[k], a_k$p
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(slot = slot_names[k], slot_p = a_k$p)
      )
    }
  }
  assert_data_frame(data, "data")
  n <- nrow(data)
  needed_vars <- character(0)
  any_W <- FALSE
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]
    if (!is.null(a_k$a)) {
      needed_vars <- c(needed_vars, all.vars(a_k$a))
    }
    if (!is.null(a_k$b)) {
      needed_vars <- c(needed_vars, all.vars(a_k$b))
    }
    if (!is.null(a_k$W)) {
      any_W <- TRUE
    }
    if (!is.null(a_k$x_vars)) {
      needed_vars <- c(needed_vars, a_k$x_vars)
    }
  }
  needed_vars <- unique(needed_vars)
  needed_vars <- intersect(needed_vars, colnames(data))
  if (length(needed_vars) > 0L) {
    has_na <- vapply(needed_vars,
                     function(v) any(is.na(data[[v]])),
                     logical(1L))
    if (any(has_na)) {
      gdpar_abort(
        sprintf(
          "Covariate(s) used by the K-individual AMM specification contain missing values: %s. Path 1 does not impute.",
          paste(sQuote(needed_vars[has_na]), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
  }
  Z_a_k_list <- vector("list", K)
  Z_b_k_list <- vector("list", K)
  Z_a_k_means_list <- vector("list", K)
  Z_b_k_means_list <- vector("list", K)
  Z_a_k_names_list <- vector("list", K)
  Z_b_k_names_list <- vector("list", K)
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]$a
    b_k <- amm_list_canonical[[k]]$b
    if (!is.null(a_k)) {
      Z_full <- stats::model.matrix(stats::update(a_k, ~ . + 0), data)
      Z_full <- as.matrix(Z_full)
      means <- colMeans(Z_full)
      Z_a_k_list[[k]] <- sweep(Z_full, 2L, means, "-")
      Z_a_k_means_list[[k]] <- means
      Z_a_k_names_list[[k]] <- colnames(Z_full)
    } else {
      Z_a_k_list[[k]] <- matrix(numeric(0), nrow = n, ncol = 0L)
      Z_a_k_means_list[[k]] <- numeric(0)
      Z_a_k_names_list[[k]] <- character(0)
    }
    if (!is.null(b_k)) {
      Z_full <- stats::model.matrix(stats::update(b_k, ~ . + 0), data)
      Z_full <- as.matrix(Z_full)
      means <- colMeans(Z_full)
      Z_b_k_list[[k]] <- sweep(Z_full, 2L, means, "-")
      Z_b_k_means_list[[k]] <- means
      Z_b_k_names_list[[k]] <- colnames(Z_full)
    } else {
      Z_b_k_list[[k]] <- matrix(numeric(0), nrow = n, ncol = 0L)
      Z_b_k_means_list[[k]] <- numeric(0)
      Z_b_k_names_list[[k]] <- character(0)
    }
  }
  if (any_W) {
    explicit_x_vars <- unique(unlist(lapply(amm_list_canonical,
                                            function(a) a$x_vars)))
    explicit_x_vars <- explicit_x_vars[nzchar(explicit_x_vars)]
    x_vars <- if (length(explicit_x_vars) > 0L) {
      explicit_x_vars
    } else {
      union_vars <- character(0)
      for (k in seq_len(K)) {
        a_k <- amm_list_canonical[[k]]
        if (!is.null(a_k$a)) {
          union_vars <- c(union_vars, all.vars(a_k$a))
        }
        if (!is.null(a_k$b)) {
          union_vars <- c(union_vars, all.vars(a_k$b))
        }
      }
      union_vars <- unique(union_vars)
      if (length(union_vars) == 0L) {
        if (inherits(formula_rhs, "formula")) {
          union_vars <- attr(stats::terms(formula_rhs, data = data),
                             "term.labels")
        } else if (is.character(formula_rhs)) {
          union_vars <- formula_rhs
        }
      }
      union_vars
    }
    if (length(x_vars) == 0L) {
      gdpar_abort(
        paste(
          "The modulating component W is active on at least one slot but no",
          "covariates were identified for the linear factor x. Include",
          "covariates in the slot's `a()`/`b()` wrappers, or supply",
          "`x_vars` to the amm_spec."
        ),
        class = "gdpar_input_error"
      )
    }
    missing_x <- setdiff(x_vars, colnames(data))
    if (length(missing_x) > 0L) {
      gdpar_abort(
        sprintf(
          "Variables required by the modulating component are missing in data: %s.",
          paste(sQuote(missing_x), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
    X_full <- as.matrix(data[, x_vars, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_means <- colMeans(X_full)
    X_centered <- sweep(X_full, 2L, X_means, "-")
    X_sds <- apply(X_centered, 2L, stats::sd)
    if (any(X_sds == 0)) {
      bad <- x_vars[X_sds == 0]
      gdpar_abort(
        sprintf(
          "Constant covariate(s) cannot enter the modulating component: %s.",
          paste(sQuote(bad), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
    X <- sweep(X_centered, 2L, X_sds, "/")
  } else {
    X <- matrix(numeric(0), nrow = n, ncol = 0L)
    X_means <- numeric(0)
    X_sds <- numeric(0)
    x_vars <- character(0)
  }
  list(
    K                = K,
    slot_names       = slot_names,
    Z_a_k_list       = Z_a_k_list,
    Z_b_k_list       = Z_b_k_list,
    X                = X,
    Z_a_k_means_list = Z_a_k_means_list,
    Z_b_k_means_list = Z_b_k_means_list,
    X_means          = X_means,
    X_sds            = X_sds,
    Z_a_k_names_list = Z_a_k_names_list,
    Z_b_k_names_list = Z_b_k_names_list,
    X_names          = if (length(x_vars) > 0L) x_vars else character(0)
  )
}

#' Print method for amm_spec objects
#'
#' @param x An object of class \code{amm_spec}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.amm_spec <- function(x, ...) {
  cat("<amm_spec> AMM Level ", x$level, "\n", sep = "")
  cat("  p (dim theta_i)    : ", x$p, "\n", sep = "")
  if (isTRUE(x$p == 1L)) {
    cat("  a (additive)       : ",
        if (is.null(x$a)) "NULL" else deparse(x$a),
        "\n", sep = "")
    cat("  b (multiplicative) : ",
        if (is.null(x$b)) "NULL" else deparse(x$b),
        "\n", sep = "")
  } else {
    cat("  dims (per-k a, b)  :\n")
    for (k in seq_len(x$p)) {
      a_k <- x$dims[[k]][["a"]]
      b_k <- x$dims[[k]][["b"]]
      cat("    k = ", k,
          " : a = ", if (is.null(a_k)) "NULL" else deparse(a_k),
          " ; b = ", if (is.null(b_k)) "NULL" else deparse(b_k),
          "\n", sep = "")
    }
  }
  cat("  W (modulating)     : ",
      if (is.null(x$W)) "NULL" else paste0("W_basis(type = '",
                                            x$W$type, "')"),
      "\n", sep = "")
  if (!is.null(x$x_vars)) {
    cat("  x_vars             : ",
        paste(x$x_vars, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

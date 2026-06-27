#' Check basis-restricted functional independence via Gram matrix
#'
#' Diagnose, before model fitting, whether the chosen finite parametric
#' representation of the AMM canonical form satisfies the
#' basis-restricted Functional Independence Condition at a candidate
#' value of the population reference. Failure indicates that the model
#' parameters are not identifiable in the chosen basis and that fitting
#' should not proceed without revising the specification.
#'
#' @param amm An object of class \code{amm_spec} produced by
#'   \code{\link{amm_spec}}, defining the bases for the additive,
#'   multiplicative and modulating components.
#' @param data A data frame containing the variables referenced in
#'   \code{amm}. Covariates are centered internally before the Gram
#'   matrix is computed, consistent with assumption (C1) of Block 1.
#' @param theta_ref_init Numeric vector of length \code{p} (the
#'   dimension of the population reference) at which the diagnostic is
#'   computed. Defaults to a vector of zeros, which corresponds to the
#'   prior mean under the default \code{\link{gdpar_prior}} on the
#'   linear-predictor scale.
#' @param formula_rhs Optional formula or character vector identifying
#'   the covariates that enter the modulating component as the linear
#'   factor x. Defaults to \code{amm$x_vars}.
#' @param tol Numeric scalar with the tolerance for the relative
#'   condition number criterion. Defaults to \code{1e-8}. The
#'   diagnostic flags failure when the smallest eigenvalue of the
#'   normalized Gram matrix is below \code{tol} times the largest
#'   eigenvalue.
#' @param family Optional \code{gdpar_family} or
#'   \code{gdpar_family_multi} object. When supplied, the report
#'   includes the parameter-level identifiability (D-ID) pre-fit layer
#'   of Block 8 Session 1 decision (3C): for each individual-scope
#'   parameter declared by the family's \code{param_specs} (Block 8
#'   Session 1 decision 1C), the declarative \code{did_status} is
#'   reported, and when the family declares K >= 2 individual
#'   parameters the symbolic separability between them is checked under
#'   \code{rigor = "full"}. In Block 8.0 K = 1 in every built-in family,
#'   so the layer reduces to a status echo; the slot is exposed for the
#'   K > 1 extension scheduled for Block 8.1. Defaults to \code{NULL}
#'   (slot omitted, backward compatible).
#' @param rigor Character scalar in \code{c("full", "fast")}
#'   controlling the C4-bis cross-component check for multivariate
#'   (\code{p > 1}) specs. Both modes (a) check the per-coordinate
#'   rank of \eqn{Z_a[k]} via its normalized Gram condition number,
#'   and (b) detect structural overlap between the column names of
#'   \eqn{Z_a[k]} and the modulating \code{X} columns (the
#'   \code{x_vars}). \code{"full"} (default) FAILS the report when
#'   any overlap is detected (conservative: structural overlap is a
#'   necessary condition for cross-component non-identifiability and
#'   the user should redesign the spec). \code{"fast"} emits a
#'   structured warning on overlap but does not flag failure (use
#'   when the user has explicitly accepted the overlap, e.g. via
#'   regularization). Note: the pre-fit check cannot detect cross-
#'   component non-identifiabilities that arise only through the
#'   posterior geometry (e.g. \eqn{\theta_{ref}}-mediated coupling
#'   between \code{a} and \code{W}); those are surfaced post-fit via
#'   divergences, low ESS and high R-hat. Ignored when \code{p == 1}.
#'
#' @return An object of class \code{gdpar_identifiability_report} with
#'   components \code{passed} (logical), \code{lambda_min},
#'   \code{lambda_max}, \code{condition_number},
#'   \code{collinear_directions} (a list describing basis-function
#'   combinations corresponding to near-zero eigenvectors when the
#'   diagnostic fails, and NULL otherwise),
#'   \code{theta_ref_used}, \code{tol_used} and \code{column_labels}.
#'   A \code{print} method provides a human-readable summary.
#'
#' @details
#' Proposition 1C of Block 1 establishes that, in a chosen finite
#' basis B for the AMM components, the basis-restricted Functional
#' Independence Condition at a value theta_ref holds if and only if
#' the population Gram matrix of the extended design matrix
#' Z_n(theta_ref) is non-singular. This function computes the empirical
#' Gram matrix from the observed data and checks its condition number.
#'
#' The diagnostic is by design local to the supplied
#' \code{theta_ref_init}: it reports the condition at one point in the
#' parameter space. The basis-restricted Functional Independence
#' Condition is a property of the basis, not of the point; in practice
#' the diagnostic at the prior mean is informative because
#' non-identifiability that is structural to the basis manifests at
#' typical reference points.
#'
#' The function does not test the abstract Functional Independence
#' Condition, which is a property of the full function classes
#' F_a, F_b, F_W. When these classes are infinite-dimensional and the
#' basis B is a finite truncation, abstract failure may occur in
#' directions outside B and would not be detected here.
#'
#' @section Methodological notes:
#' The threshold criterion is relative
#' (smallest eigenvalue divided by the largest, compared against
#' \code{tol}), not absolute. The relative criterion coincides with
#' the inverse of the condition number, which is invariant to
#' rescaling of the basis columns. An absolute eigenvalue threshold
#' would depend on the scale of the covariates and would produce false
#' positives when the basis matrix has columns of very different
#' magnitudes. To make the diagnostic robust to scale heterogeneity
#' across basis terms, the columns of the extended design matrix are
#' normalized to unit norm before the Gram matrix is computed.
#' Normalization affects only the diagnostic; model fitting uses the
#' matrices in their natural scale.
#'
#' When the diagnostic fails, \code{collinear_directions} is populated
#' by projecting the eigenvectors associated with eigenvalues below
#' \code{tol} times the largest eigenvalue onto the named basis
#' columns, producing a human-readable description of which
#' combinations of basis functions are linearly dependent.
#'
#' @section Dependencies:
#' This function calls \code{\link[base]{eigen}} with
#' \code{symmetric = TRUE} for the eigendecomposition of the Gram
#' matrix, and \code{\link[stats]{model.matrix}} for evaluating the
#' formula-based bases of the additive and multiplicative components.
#'
#' @references
#' Block 1 of the package theoretical addendum, Section 6.6,
#' Proposition 1C. See
#' \code{vignette("01_amm_identifiability", package = "gdpar")}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
#' spec <- amm_spec(
#'   a = ~ x1 + x2,
#'   b = ~ x1,
#'   W = W_basis(type = "polynomial", degree = 1),
#'   x_vars = "x2"
#' )
#' report <- gdpar_check_identifiability(spec, df, theta_ref_init = 0.5)
#' print(report)
#'
#' @seealso \code{\link{amm_spec}}, \code{\link{W_basis}},
#'   \code{\link{gdpar}}
#' @export
gdpar_check_identifiability <- function(amm, data, theta_ref_init = NULL,
                                        formula_rhs = NULL,
                                        family = NULL, tol = 1e-8,
                                        rigor = c("full", "fast")) {
  assert_inherits(amm, "amm_spec", "amm")
  assert_data_frame(data, "data")
  assert_numeric_scalar(tol, "tol", lower = 0, upper = 1)
  rigor <- match.arg(rigor)
  if (amm$level == 0L) {
    obj <- list(
      passed = TRUE,
      lambda_min = NA_real_,
      lambda_max = NA_real_,
      condition_number = NA_real_,
      collinear_directions = NULL,
      theta_ref_used = theta_ref_init,
      tol_used = tol,
      column_labels = character(0),
      message = "AMM Level 0: identifiability is trivial."
    )
    class(obj) <- c("gdpar_identifiability_report", "list")
    return(obj)
  }
  if (is.null(theta_ref_init)) {
    theta_ref_init <- if (!is.null(amm$b)) 1 else 0
  }
  if (!is.numeric(theta_ref_init)) {
    gdpar_abort(
      "Argument 'theta_ref_init' must be numeric.",
      class = "gdpar_input_error"
    )
  }
  if (!is.null(amm$b) && all(theta_ref_init == 0)) {
    gdpar_warn(
      paste(
        "theta_ref_init is zero in every coordinate while the",
        "multiplicative component b is active; the diagnostic block",
        "for b vanishes trivially and the report cannot detect b-related",
        "non-identifiability. Pass a non-zero theta_ref_init for a",
        "meaningful diagnostic of the b block."
      ),
      class = "gdpar_diagnostic_warning"
    )
  }
  p <- length(theta_ref_init)
  if (p < 1L) {
    gdpar_abort(
      "Argument 'theta_ref_init' must have length at least one.",
      class = "gdpar_input_error"
    )
  }
  if (!is.null(amm$W)) {
    amm$W <- materialize_W_basis(amm$W, p)
  }
  design <- build_amm_design(amm, data, formula_rhs)
  if (is.null(design$Z_a) && !is.null(design$Z_a_list)) {
    return(check_identifiability_multi(
      amm           = amm,
      design        = design,
      theta_ref     = theta_ref_init,
      tol           = tol,
      rigor         = rigor,
      family        = family
    ))
  }
  Z_a <- design$Z_a
  Z_b <- design$Z_b
  X <- design$X
  blocks <- list()
  labels <- character(0)
  if (ncol(Z_a) > 0L) {
    blocks$a <- Z_a
    labels <- c(labels, paste0("a:", design$Z_a_names))
  }
  if (ncol(Z_b) > 0L) {
    if (length(theta_ref_init) == 1L) {
      Z_b_anchored <- Z_b * theta_ref_init
      blocks$b <- Z_b_anchored
      labels <- c(labels, paste0("b*theta:", design$Z_b_names))
    } else {
      Z_b_anchored <- do.call(
        cbind,
        lapply(seq_along(theta_ref_init), function(k) {
          Z_b * theta_ref_init[k]
        })
      )
      labs <- as.vector(outer(design$Z_b_names,
                              seq_along(theta_ref_init),
                              function(a, k) paste0("b*theta[", k, "]:", a)))
      blocks$b <- Z_b_anchored
      labels <- c(labels, labs)
    }
  }
  if (ncol(X) > 0L) {
    blocks$W <- X
    labels <- c(labels, paste0("W*x:", design$X_names))
  }
  if (length(blocks) == 0L) {
    obj <- list(
      passed = TRUE, lambda_min = NA_real_, lambda_max = NA_real_,
      condition_number = NA_real_,
      collinear_directions = NULL,
      theta_ref_used = theta_ref_init, tol_used = tol,
      column_labels = character(0),
      message = "No active blocks; identifiability check is trivial."
    )
    class(obj) <- c("gdpar_identifiability_report", "list")
    return(obj)
  }
  Z <- do.call(cbind, blocks)
  col_norms <- sqrt(colSums(Z * Z))
  zero_norm_cols <- which(col_norms == 0)
  if (length(zero_norm_cols) > 0L) {
    obj <- list(
      passed = FALSE, lambda_min = 0, lambda_max = NA_real_,
      condition_number = Inf,
      collinear_directions = list(
        list(label = "zero-norm columns",
             columns = labels[zero_norm_cols],
             coefficients = rep(1, length(zero_norm_cols)))
      ),
      theta_ref_used = theta_ref_init, tol_used = tol,
      column_labels = labels,
      message = paste(
        "Some basis columns are zero after centering.",
        "Identifiability fails in those directions."
      )
    )
    class(obj) <- c("gdpar_identifiability_report", "list")
    return(obj)
  }
  Z_normalized <- sweep(Z, 2L, col_norms, "/")
  G <- crossprod(Z_normalized) / nrow(Z_normalized)
  eig <- eigen(G, symmetric = TRUE)
  lambda_max <- max(eig$values)
  lambda_min <- min(eig$values)
  condition_number <- lambda_max / max(lambda_min, .Machine$double.eps)
  passed <- lambda_min >= tol * lambda_max
  collinear_directions <- NULL
  if (!passed) {
    bad_idx <- which(eig$values < tol * lambda_max)
    collinear_directions <- lapply(bad_idx, function(j) {
      v <- eig$vectors[, j]
      keep <- which(abs(v) > 1e-3)
      ord <- order(abs(v[keep]), decreasing = TRUE)
      list(
        eigenvalue = eig$values[j],
        columns = labels[keep][ord],
        coefficients = v[keep][ord]
      )
    })
  }
  c4_bis <- NULL
  passed_c4_bis <- TRUE
  if (p > 1L && !is.null(amm$W) && ncol(X) > 0L &&
      !is.null(design$Z_a_list)) {
    c4_bis <- check_C4_bis_per_k(
      design       = design,
      amm          = amm,
      theta_ref    = theta_ref_init,
      rigor        = rigor,
      tol          = tol
    )
    passed_c4_bis <- all(vapply(c4_bis$per_k,
                                function(z) isTRUE(z$passed),
                                logical(1L)))
  }

  did_pre_fit <- .check_did_pre_fit(family, design, theta_ref_init, rigor)
  passed_did <- if (is.null(did_pre_fit)) TRUE else isTRUE(did_pre_fit$passed)
  passed_final <- passed && passed_c4_bis && passed_did
  obj <- list(
    passed = passed_final,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    condition_number = condition_number,
    collinear_directions = collinear_directions,
    theta_ref_used = theta_ref_init,
    tol_used = tol,
    rigor_used = rigor,
    column_labels = labels,
    c4_bis = c4_bis,
    did_pre_fit = did_pre_fit,
    message = if (passed_final) {
      if (p > 1L && !is.null(c4_bis)) {
        sprintf(
          "Basis-restricted Functional Independence Condition (C1-C4 + C4-bis, rigor=%s) holds at theta_ref_init.",
          rigor
        )
      } else {
        "Basis-restricted Functional Independence Condition holds at theta_ref_init."
      }
    } else if (!passed) {
      "Basis-restricted Functional Independence Condition FAILS at theta_ref_init."
    } else if (!passed_c4_bis) {
      sprintf(
        "C1-C4 hold but C4-bis (cross-component per-k, rigor=%s) FAILS for at least one coordinate. Inspect 'c4_bis$per_k'.",
        rigor
      )
    } else {
      "C1-C4 (and C4-bis when applicable) hold, but pre-fit parameter identifiability (D-ID) FAILS for at least one individual parameter. Inspect 'did_pre_fit$per_param'."
    }
  )
  class(obj) <- c("gdpar_identifiability_report", "list")
  obj
}


#' Internal: identifiability check for the multivariate (p > 1) path
#'
#' Counterpart of the univariate Gram-matrix C1-C4 check for the
#' multivariate spec. Iterates over coordinates: for each \eqn{k}
#' builds the per-coordinate extended Gram matrix (additive design
#' \eqn{Z_a[k]} plus, when the modulating component is present, the
#' columns \eqn{(\theta_{ref}[k]^m - 0^m) X} for \eqn{m = 1, \ldots,
#' W_{per_k_dim}}), normalizes the columns and tests the condition
#' number. The aggregated report exposes \code{lambda_min},
#' \code{lambda_max}, \code{condition_number} as the worst-case
#' (smallest min eigenvalue / largest condition number) across
#' coordinates, plus the full per-k breakdown in \code{c4_bis}.
#'
#' @keywords internal
#' @noRd
check_identifiability_multi <- function(amm, design, theta_ref,
                                        tol, rigor, family = NULL) {
  p <- length(theta_ref)
  c4_bis <- check_C4_bis_per_k(
    design    = design,
    amm       = amm,
    theta_ref = theta_ref,
    rigor     = rigor,
    tol       = tol
  )
  passed_c4_bis <- all(vapply(c4_bis$per_k,
                              function(z) isTRUE(z$passed),
                              logical(1L)))
  lmins <- vapply(c4_bis$per_k,
                  function(z) z$lambda_min %||% NA_real_,
                  numeric(1L))
  lmaxs <- vapply(c4_bis$per_k,
                  function(z) z$lambda_max %||% NA_real_,
                  numeric(1L))
  conds <- vapply(c4_bis$per_k,
                  function(z) z$condition_number %||% NA_real_,
                  numeric(1L))
  did_pre_fit <- .check_did_pre_fit(family, design, theta_ref, rigor)
  passed_did <- if (is.null(did_pre_fit)) TRUE else isTRUE(did_pre_fit$passed)
  passed_all <- passed_c4_bis && passed_did
  obj <- list(
    passed               = passed_all,
    lambda_min           = if (all(is.na(lmins))) NA_real_
                           else min(lmins, na.rm = TRUE),
    lambda_max           = if (all(is.na(lmaxs))) NA_real_
                           else max(lmaxs, na.rm = TRUE),
    condition_number     = if (all(is.na(conds))) NA_real_
                           else max(conds, na.rm = TRUE),
    collinear_directions = NULL,
    theta_ref_used       = theta_ref,
    tol_used             = tol,
    rigor_used           = rigor,
    column_labels        = character(0),
    c4_bis               = c4_bis,
    did_pre_fit          = did_pre_fit,
    message              = if (passed_all) {
      sprintf(
        "Multivariate basis-restricted Functional Independence Condition (C4-bis per-k, rigor=%s) holds at theta_ref_init across all p=%d coordinates.",
        rigor, p
      )
    } else if (!passed_c4_bis) {
      sprintf(
        "Multivariate basis-restricted Functional Independence Condition (C4-bis per-k, rigor=%s) FAILS for at least one coordinate. Inspect 'c4_bis$per_k'.",
        rigor
      )
    } else {
      sprintf(
        "Pre-fit parameter identifiability (D-ID) layer FAILS for at least one individual parameter. Inspect 'did_pre_fit$per_param'."
      )
    }
  )
  class(obj) <- c("gdpar_identifiability_report", "list")
  obj
}


#' Internal: C4-bis per-coordinate cross-component identifiability check
#'
#' For each coordinate \eqn{k = 1, \ldots, p}, builds the extended
#' Gram matrix from the additive design \eqn{Z_a[k]} and the
#' modulating-component design \eqn{(\theta_{ref}[k]^m -
#' \theta_{anchor}[k]^m) X} for \eqn{m = 1, \ldots, W_{per_k_dim}}
#' (the anchor is taken as zero, consistent with the package default
#' on the linear-predictor scale), normalizes the columns to unit
#' norm and computes the condition number of the normalized Gram.
#' \code{rigor = "fast"} skips the Gram construction and only checks
#' the structural overlap between the column names of \eqn{Z_a[k]}
#' and \eqn{X}; when an overlap is detected the function emits a
#' \code{gdpar_c4bis_overlap_warning} but does not flag failure.
#'
#' @param design List returned by \code{build_amm_design}; expected to
#'   contain \code{Z_a_list} (per-coordinate additive design),
#'   \code{Z_a_names_list} (per-coordinate column labels), \code{X}
#'   (modulating-component design), and \code{X_names}.
#' @param amm The \code{amm_spec}; used to read \code{W$dim} and
#'   \code{p} to compute \code{W_per_k_dim}.
#' @param theta_ref Numeric vector of length \code{p}.
#' @param rigor Character scalar (\code{"full"} or \code{"fast"}).
#' @param tol Numeric tolerance for the condition number criterion.
#' @return A list with elements \code{rigor} (echoed) and \code{per_k}
#'   (a list of length \code{p} with per-coordinate diagnostics).
#'   Each per-k entry has \code{passed}, \code{lambda_min},
#'   \code{lambda_max}, \code{condition_number}, \code{shared_cols}
#'   (character vector of overlapping column names),
#'   \code{collinear_directions}.
#' @keywords internal
#' @noRd
check_C4_bis_per_k <- function(design, amm, theta_ref, rigor, tol) {
  p <- length(theta_ref)
  per_k_dim <- if (!is.null(amm$W$dim) && p >= 1L) {
    as.integer(amm$W$dim / p)
  } else {
    0L
  }
  X <- design$X
  X_names <- design$X_names %||% colnames(X) %||% character(ncol(X))
  Z_a_list <- design$Z_a_list
  Z_a_names_list <- design$Z_a_names_list
  per_k <- vector("list", p)
  overlap_acc <- vector("list", p)  # accumulate for consolidated warning
  for (k in seq_len(p)) {
    Z_a_k <- if (length(Z_a_list) >= k) Z_a_list[[k]] else matrix(0, 0, 0)
    Z_a_k_names <- if (length(Z_a_names_list) >= k) {
      Z_a_names_list[[k]]
    } else if (!is.null(colnames(Z_a_k))) {
      colnames(Z_a_k)
    } else {
      character(ncol(Z_a_k))
    }
    shared <- intersect(Z_a_k_names, X_names)
    if (length(shared) > 0L) {
      overlap_acc[[k]] <- shared
    }

    if (identical(rigor, "fast")) {
      per_k[[k]] <- list(
        passed               = TRUE,
        rigor                = "fast",
        lambda_min           = NA_real_,
        lambda_max           = NA_real_,
        condition_number     = NA_real_,
        shared_cols          = shared,
        collinear_directions = NULL,
        coord                = k
      )
      next
    }

    # rigor == "full": rank check of Z_a[k] alone (per-k version of
    # C1-C4) PLUS structural overlap check between Z_a[k] columns
    # and X columns. Rationale: building the extended Gram matrix
    # [Z_a[k] | (theta_ref[k]^m) X]_{m=1..W_per_k_dim} produces
    # spurious rank deficit because at a fixed theta_ref the W
    # columns are all scalar multiples of X per coordinate l, so
    # rank(W-block per l) = 1 by construction. A Jacobian-based
    # check at (a = 0, W = 0) inherits the same defect because the
    # information of theta_ref in W*x evaporates when W = 0. The
    # rigorous detection of cross-component non-identifiability
    # requires the posterior geometry (divergences, low ESS, high
    # R-hat) and lives in the post-fit diagnostics, not in this
    # pre-fit structural check.
    if (ncol(Z_a_k) == 0L) {
      per_k[[k]] <- list(
        passed               = TRUE,
        rigor                = "full",
        lambda_min           = NA_real_,
        lambda_max           = NA_real_,
        condition_number     = NA_real_,
        shared_cols          = shared,
        collinear_directions = NULL,
        coord                = k,
        message              = "empty Z_a[k]: trivially passes."
      )
      next
    }
    col_norms_k <- sqrt(colSums(Z_a_k * Z_a_k))
    if (any(col_norms_k == 0)) {
      bad <- which(col_norms_k == 0)
      per_k[[k]] <- list(
        passed               = FALSE,
        rigor                = "full",
        lambda_min           = 0,
        lambda_max           = NA_real_,
        condition_number     = Inf,
        shared_cols          = shared,
        collinear_directions = list(
          list(label        = "zero-norm columns in Z_a[k]",
               columns      = paste0("a", k, ":",
                                     Z_a_k_names[bad]),
               coefficients = rep(1, length(bad)))
        ),
        coord                = k
      )
      next
    }
    Zk_n <- sweep(Z_a_k, 2L, col_norms_k, "/")
    Gk <- crossprod(Zk_n) / nrow(Zk_n)
    eig_k <- eigen(Gk, symmetric = TRUE)
    lmax <- max(eig_k$values)
    lmin <- min(eig_k$values)
    cond_k <- lmax / max(lmin, .Machine$double.eps)
    passed_rank <- lmin >= tol * lmax
    passed_overlap <- length(shared) == 0L
    passed_k <- passed_rank && passed_overlap
    cdirs <- NULL
    if (!passed_rank) {
      bad_idx <- which(eig_k$values < tol * lmax)
      cdirs <- lapply(bad_idx, function(j) {
        v <- eig_k$vectors[, j]
        keep <- which(abs(v) > 1e-3)
        ord <- order(abs(v[keep]), decreasing = TRUE)
        list(
          eigenvalue   = eig_k$values[j],
          columns      = paste0("a", k, ":",
                                Z_a_k_names[keep][ord]),
          coefficients = v[keep][ord]
        )
      })
    }
    if (!passed_overlap) {
      cdirs <- c(
        cdirs,
        list(list(
          label        = sprintf(
            "C4-bis structural overlap on coord %d", k),
          columns      = paste0("shared with X: ", shared),
          coefficients = rep(NA_real_, length(shared))
        ))
      )
    }
    per_k[[k]] <- list(
      passed               = passed_k,
      rigor                = "full",
      lambda_min           = lmin,
      lambda_max           = lmax,
      condition_number     = cond_k,
      shared_cols          = shared,
      collinear_directions = cdirs,
      coord                = k
    )
  }
  if (identical(rigor, "fast")) {
    overlap_idx <- which(!vapply(overlap_acc, is.null, logical(1L)))
    if (length(overlap_idx) > 0L) {
      detail <- vapply(overlap_idx, function(k) {
        sprintf("coord %d: %s",
                k,
                paste(sQuote(overlap_acc[[k]]), collapse = ", "))
      }, character(1L))
      gdpar_warn(
        paste0(
          "C4-bis (rigor = 'fast'): structural overlap between ",
          "additive component a and modulating x_vars detected for ",
          length(overlap_idx),
          " coordinate(s) [",
          paste(detail, collapse = "; "),
          "]. Overlap is necessary but not sufficient for ",
          "non-identifiability; run with rigor = 'full' to apply ",
          "the conservative check."
        ),
        class = "gdpar_c4bis_overlap_warning",
        data = list(
          shared_cols_by_coord = setNames(
            overlap_acc[overlap_idx],
            paste0("coord_", overlap_idx)
          )
        )
      )
    }
  }
  list(rigor = rigor, per_k = per_k)
}

#' Internal: pre-fit parameter identifiability (D-ID) layer
#'
#' Block 8 Session 1 decision (3C), layer (ii): for each individual-scope
#' parameter (those with \code{scope} equal to \code{"per_observation"}
#' or \code{"per_group"}) of the family's \code{param_specs}, report the
#' declarative \code{did_status}, and for K >= 2 check symbolic
#' separability between parameters under \code{rigor = "full"}.
#'
#' In Block 8.0 every built-in family declares K = 1 (only the first
#' \code{param_spec} is individual-scope; remaining entries are
#' \code{scope = "population"} auxiliaries), so the layer reduces to
#' echoing the K = 1 status. The K > 1 extension queued for Block 8.1
#' will exercise the symbolic-separability sub-check; the slot
#' \code{symbolic_separability} is exposed already but populated only
#' when K >= 2.
#'
#' Symbolic separability is interpreted minimally in Block 8.0: two
#' individual params with the same \code{prior_canonical_kind} are
#' flagged as non-separable (their roles overlap canonically) unless
#' the user has explicitly distinguished them via different
#' \code{name}. This check is conservative; finer separability analyses
#' are queued for the K > 1 implementation block.
#'
#' @param family Either NULL (slot omitted, function returns NULL),
#'   a \code{gdpar_family}, or a \code{gdpar_family_multi}. In the
#'   multivariate case the family is expected to be homogeneous across
#'   coordinates; the per-parameter layer is evaluated on
#'   \code{families[[1L]]} and applies to every coordinate identically.
#' @param design Ignored in Block 8.0 (slot for the rank-by-parameter
#'   sub-check queued for Block 8.1, which will read \code{Z_a_list} /
#'   \code{Z_b_list} per individual parameter).
#' @param theta_ref Numeric vector with the reference point used in the
#'   surrounding C1-C4 check.
#' @param rigor Character scalar; only \code{"full"} runs the symbolic
#'   separability sub-check (no-op for K = 1).
#' @return Either NULL (when family is NULL) or a list with components
#'   \code{passed} (logical), \code{K} (integer, count of individual
#'   parameters), \code{per_param} (list of length K with
#'   \code{name}, \code{scope}, \code{did_status}, \code{did_condition},
#'   \code{did_reference}, \code{prior_canonical_kind}), and
#'   \code{symbolic_separability} (list; populated only when K >= 2).
#' @keywords internal
#' @noRd
.check_did_pre_fit <- function(family, design, theta_ref, rigor) {
  if (is.null(family)) {
    return(NULL)
  }
  if (inherits(family, "gdpar_family_multi")) {
    base_family <- family$families[[1L]]
  } else if (inherits(family, "gdpar_family")) {
    base_family <- family
  } else {
    return(NULL)
  }
  if (is.null(base_family$param_specs)) {
    return(NULL)
  }
  individual_scopes <- c("per_observation", "per_group")
  is_individual <- vapply(base_family$param_specs,
                          function(s) s$scope %in% individual_scopes,
                          logical(1L))
  individual_specs <- base_family$param_specs[is_individual]
  K <- length(individual_specs)
  per_param <- lapply(individual_specs, function(s) {
    list(
      name                 = s$name,
      scope                = s$scope,
      did_status           = s$did_status,
      did_condition        = s$did_condition,
      did_reference        = s$did_reference,
      prior_canonical_kind = s$prior_canonical_kind
    )
  })
  passed_did_declarative <- all(vapply(individual_specs,
    function(s) s$did_status %in% c("holds", "holds_under_condition"),
    logical(1L)
  ))
  symbolic_separability <- NULL
  passed_separability <- TRUE
  if (K >= 2L && identical(rigor, "full")) {
    prior_kinds <- vapply(individual_specs,
                          function(s) s$prior_canonical_kind,
                          character(1L))
    names_per <- vapply(individual_specs, function(s) s$name,
                        character(1L))
    overlap <- duplicated(prior_kinds) | duplicated(prior_kinds,
                                                    fromLast = TRUE)
    if (any(overlap)) {
      offending <- which(overlap)
      symbolic_separability <- list(
        passed = FALSE,
        overlapping_kinds = prior_kinds[offending],
        overlapping_names = names_per[offending],
        message = paste(
          "Symbolic separability: two or more individual parameters",
          "declare the same prior_canonical_kind. Their identification",
          "from data may collapse without an explicit canonical",
          "distinction between them."
        )
      )
      passed_separability <- FALSE
    } else {
      symbolic_separability <- list(
        passed = TRUE,
        overlapping_kinds = character(0),
        overlapping_names = character(0),
        message = "Symbolic separability: distinct prior_canonical_kind per parameter."
      )
    }
  }
  list(
    passed                = passed_did_declarative && passed_separability,
    K                     = K,
    per_param             = per_param,
    symbolic_separability = symbolic_separability,
    rigor                 = rigor
  )
}

#' Internal: per-slot rank check of Z_a per slot for the K-individual path
#'
#' D-B3 of sub-phase 8.3.4 (Block 8): pre-fit structural rank check of
#' the per-slot additive design matrix Z_a per slot for each
#' K-individual slot k = 1, ..., K. Each Z_a per slot is
#' column-centered by \code{.build_amm_design_K()} before reaching
#' this helper; the check computes the normalized Gram condition
#' number per slot and flags slots whose minimum eigenvalue falls
#' below \code{tol * lambda_max} as non-identifiable.
#'
#' Complementary to D-B2 (cross-slot extended Gram check) and the
#' per-slot univariate \code{gdpar_check_identifiability()} that is
#' wired upstream in \code{.gdpar_K()}. The three layers together
#' cover (i) per-slot AMM canonical (C1-C4), (ii) per-slot Z_a rank
#' (this helper, D-B3), (iii) cross-slot rank (.check_C4_bis_K, D-B2).
#'
#' @param design_K List returned by \code{.build_amm_design_K()} with
#'   components \code{Z_a_k_list}, \code{Z_a_k_names_list},
#'   \code{slot_names}.
#' @param rigor Character scalar; \code{"fast"} returns NA diagnostics
#'   (rank is checked only structurally via column-zero detection);
#'   \code{"full"} computes the Gram eigendecomposition.
#' @param tol Numeric tolerance for the condition criterion. The slot
#'   passes when \code{lambda_min >= tol * lambda_max}.
#' @return List with \code{passed} (logical scalar), \code{rigor},
#'   \code{per_slot} (list of length K).
#' @keywords internal
#' @noRd
.check_Z_a_K_per_slot <- function(design_K, rigor, tol) {
  Z_a_list <- design_K$Z_a_k_list
  names_list <- design_K$Z_a_k_names_list
  slot_names <- design_K$slot_names
  K <- length(Z_a_list)
  per_slot <- vector("list", K)
  passed_all <- TRUE
  for (k in seq_len(K)) {
    Z_a_k <- Z_a_list[[k]]
    z_names <- if (length(names_list) >= k && !is.null(names_list[[k]])) {
      names_list[[k]]
    } else if (!is.null(colnames(Z_a_k))) {
      colnames(Z_a_k)
    } else {
      character(ncol(Z_a_k))
    }
    if (ncol(Z_a_k) == 0L) {
      per_slot[[k]] <- list(
        slot              = slot_names[k],
        passed            = TRUE,
        rigor             = rigor,
        lambda_min        = NA_real_,
        lambda_max        = NA_real_,
        condition_number  = NA_real_,
        collinear_columns = NULL,
        message           = "empty Z_a^{(k)}: trivially passes"
      )
      next
    }
    col_norms <- sqrt(colSums(Z_a_k * Z_a_k))
    if (any(col_norms == 0)) {
      bad <- which(col_norms == 0)
      per_slot[[k]] <- list(
        slot              = slot_names[k],
        passed            = FALSE,
        rigor             = rigor,
        lambda_min        = 0,
        lambda_max        = NA_real_,
        condition_number  = Inf,
        collinear_columns = z_names[bad],
        message           = sprintf(
          "Z_a^{(%s)} has %d zero-norm column(s): %s",
          slot_names[k], length(bad),
          paste(sQuote(z_names[bad]), collapse = ", ")
        )
      )
      passed_all <- FALSE
      next
    }
    if (identical(rigor, "fast")) {
      per_slot[[k]] <- list(
        slot              = slot_names[k],
        passed            = TRUE,
        rigor             = "fast",
        lambda_min        = NA_real_,
        lambda_max        = NA_real_,
        condition_number  = NA_real_,
        collinear_columns = NULL,
        message           = "rigor='fast' skipped Gram eigendecomposition"
      )
      next
    }
    Zn <- sweep(Z_a_k, 2L, col_norms, "/")
    G <- crossprod(Zn) / nrow(Zn)
    eg <- eigen(G, symmetric = TRUE)
    lmin <- min(eg$values)
    lmax <- max(eg$values)
    cnd <- lmax / max(lmin, .Machine$double.eps)
    ok <- lmin >= tol * lmax
    cdirs <- NULL
    if (!ok) {
      bad_idx <- which(eg$values < tol * lmax)
      cdirs <- lapply(bad_idx, function(j) {
        v <- eg$vectors[, j]
        keep <- which(abs(v) > 1e-3)
        ord <- order(abs(v[keep]), decreasing = TRUE)
        list(
          eigenvalue   = eg$values[j],
          columns      = z_names[keep][ord],
          coefficients = v[keep][ord]
        )
      })
      passed_all <- FALSE
    }
    per_slot[[k]] <- list(
      slot              = slot_names[k],
      passed            = ok,
      rigor             = "full",
      lambda_min        = lmin,
      lambda_max        = lmax,
      condition_number  = cnd,
      collinear_columns = cdirs,
      message           = if (ok) {
        "Z_a^{(k)} full column rank under tol criterion"
      } else {
        sprintf(
          paste(
            "Z_a^{(%s)} rank-deficient: lambda_min/lambda_max = %.3e",
            "below tol = %.3e"
          ),
          slot_names[k], lmin / max(lmax, .Machine$double.eps), tol
        )
      }
    )
  }
  names(per_slot) <- slot_names
  list(passed = passed_all, rigor = rigor, per_slot = per_slot)
}

#' Internal: cross-slot extended Gram check for the K-individual path
#'
#' D-B2 of sub-phase 8.3.4 (Block 8): pre-fit structural rank check on
#' the column-wise concatenation of the per-slot additive design
#' matrices Z_a^(1), ..., Z_a^(K). Even when each per-slot matrix alone
#' is full column
#' rank (D-B3), the joint matrix can be rank-deficient when the same
#' covariate appears in multiple slots with linearly equivalent
#' designs. This check flags cross-slot collinearity that escapes the
#' per-slot D-B3 layer.
#'
#' Naming: each column of the concatenation is prefixed with the slot
#' name, so collinear-direction reports identify which slot
#' contributes each column.
#'
#' @param design_K List returned by \code{.build_amm_design_K()}.
#' @param rigor Character scalar (\code{"fast"} or \code{"full"});
#'   \code{"fast"} returns a structural pass without eigendecomposition.
#' @param tol Numeric tolerance for the condition criterion.
#' @return List with \code{passed} (logical), \code{rigor},
#'   \code{lambda_min}, \code{lambda_max}, \code{condition_number},
#'   \code{collinear_directions}, \code{total_columns}.
#' @keywords internal
#' @noRd
.check_C4_bis_K_cross_slot <- function(design_K, rigor, tol) {
  Z_a_list <- design_K$Z_a_k_list
  names_list <- design_K$Z_a_k_names_list
  slot_names <- design_K$slot_names
  K <- length(Z_a_list)
  prefixed_cols <- list()
  combined_blocks <- list()
  for (k in seq_len(K)) {
    Z_a_k <- Z_a_list[[k]]
    if (ncol(Z_a_k) == 0L) next
    z_names <- if (length(names_list) >= k && !is.null(names_list[[k]])) {
      names_list[[k]]
    } else if (!is.null(colnames(Z_a_k))) {
      colnames(Z_a_k)
    } else {
      paste0("col", seq_len(ncol(Z_a_k)))
    }
    prefixed_cols[[length(prefixed_cols) + 1L]] <-
      paste0(slot_names[k], ":", z_names)
    combined_blocks[[length(combined_blocks) + 1L]] <- Z_a_k
  }
  if (length(combined_blocks) == 0L) {
    return(list(
      passed               = TRUE,
      rigor                = rigor,
      lambda_min           = NA_real_,
      lambda_max           = NA_real_,
      condition_number     = NA_real_,
      collinear_directions = NULL,
      total_columns        = 0L,
      message              = "empty cross-slot design"
    ))
  }
  Z_joint <- do.call(cbind, combined_blocks)
  joint_names <- unlist(prefixed_cols)
  if (identical(rigor, "fast")) {
    return(list(
      passed               = TRUE,
      rigor                = "fast",
      lambda_min           = NA_real_,
      lambda_max           = NA_real_,
      condition_number     = NA_real_,
      collinear_directions = NULL,
      total_columns        = ncol(Z_joint),
      message              = "rigor='fast' skipped cross-slot eigendecomposition"
    ))
  }
  col_norms <- sqrt(colSums(Z_joint * Z_joint))
  if (any(col_norms == 0)) {
    bad <- which(col_norms == 0)
    return(list(
      passed               = FALSE,
      rigor                = "full",
      lambda_min           = 0,
      lambda_max           = NA_real_,
      condition_number     = Inf,
      collinear_directions = list(list(
        eigenvalue   = 0,
        columns      = joint_names[bad],
        coefficients = rep(1, length(bad))
      )),
      total_columns        = ncol(Z_joint),
      message              = "zero-norm columns in cross-slot design"
    ))
  }
  Zn <- sweep(Z_joint, 2L, col_norms, "/")
  G <- crossprod(Zn) / nrow(Zn)
  eg <- eigen(G, symmetric = TRUE)
  lmin <- min(eg$values)
  lmax <- max(eg$values)
  cnd <- lmax / max(lmin, .Machine$double.eps)
  ok <- lmin >= tol * lmax
  cdirs <- NULL
  if (!ok) {
    bad_idx <- which(eg$values < tol * lmax)
    cdirs <- lapply(bad_idx, function(j) {
      v <- eg$vectors[, j]
      keep <- which(abs(v) > 1e-3)
      ord <- order(abs(v[keep]), decreasing = TRUE)
      list(
        eigenvalue   = eg$values[j],
        columns      = joint_names[keep][ord],
        coefficients = v[keep][ord]
      )
    })
  }
  list(
    passed               = ok,
    rigor                = "full",
    lambda_min           = lmin,
    lambda_max           = lmax,
    condition_number     = cnd,
    collinear_directions = cdirs,
    total_columns        = ncol(Z_joint),
    message              = if (ok) {
      "cross-slot extended Gram full rank under tol criterion"
    } else {
      sprintf(
        "cross-slot extended Gram rank-deficient (cond = %.3e, tol = %.3e)",
        cnd, tol
      )
    }
  )
}

#' Internal: K-individual identifiability check (D-B2 + D-B3)
#'
#' Combines the per-slot rank check (D-B3) and the cross-slot extended
#' Gram check (D-B2) into a single structured report consumed by
#' \code{.gdpar_K()} before the Stan fit. Layered with the per-slot
#' univariate \code{gdpar_check_identifiability()} that
#' \code{.gdpar_K()} also runs (which exercises C1-C4 within each
#' slot's univariate AMM).
#'
#' Sub-phase 8.3.4 of Block 8 introduces this entry. The K = 1 path
#' continues to use \code{gdpar_check_identifiability()} directly; this
#' helper is only called when K >= 2.
#'
#' @param design_K List returned by \code{.build_amm_design_K()}.
#' @param rigor Character scalar; \code{"fast"} or \code{"full"}.
#' @param tol Numeric tolerance for the rank criteria.
#' @return List with \code{passed} (logical, conjunction over layers),
#'   \code{rigor}, \code{tol}, \code{per_slot_rank} (D-B3 report),
#'   \code{cross_slot_gram} (D-B2 report), \code{K}.
#' @keywords internal
#' @noRd
.check_identifiability_K <- function(design_K, rigor = "full", tol = 1e-8) {
  rigor <- match.arg(rigor, c("full", "fast"))
  per_slot_rank <- .check_Z_a_K_per_slot(design_K, rigor = rigor,
                                          tol = tol)
  cross_slot_gram <- .check_C4_bis_K_cross_slot(design_K, rigor = rigor,
                                                 tol = tol)
  list(
    passed          = per_slot_rank$passed && cross_slot_gram$passed,
    rigor           = rigor,
    tol             = tol,
    K               = length(design_K$Z_a_k_list),
    per_slot_rank   = per_slot_rank,
    cross_slot_gram = cross_slot_gram
  )
}

#' Internal: post-fit information contraction per K-individual slot
#'
#' D-B1 of sub-phase 8.3.4 (Block 8): for each per-slot anchor
#' \code{theta_ref_k[, k]} (single-anchor regime) or
#' \code{mu_theta_ref_k[1, k]} (per-group regime), compute the
#' prior-to-posterior contraction
#' \code{C_k = 1 - var_post(theta_ref_k) / var_prior(theta_ref_k)}.
#' The prior variance is read from the canonical kind associated with
#' the slot's \code{prior_canonical_kind} via
#' \code{.gdpar_canonical_prior_variance()}; if the user has overridden
#' the prior via \code{priors_by_kind} the override is silently used as
#' long as it is parseable as a Stan \code{normal(loc, scale)} or
#' \code{student_t(df, loc, scale)} expression (Block 8.3.4 supports
#' these two; other kinds fall back to the canonical default with a
#' diagnostic message).
#'
#' Thresholds (decision accepted in the sub-phase 8.3.4 scoping):
#' \itemize{
#'   \item \code{C_k >= 0.5}: pass (data has substantially informed
#'         the slot anchor).
#'   \item \code{0.1 <= C_k < 0.5}: warn (weak learning; the prior
#'         still dominates the posterior by more than half its
#'         variance).
#'   \item \code{C_k < 0.1}: information-error (data is essentially
#'         uninformative for this slot anchor; the posterior recovers
#'         the prior almost exactly).
#' }
#' The information-error is reported via the diagnostics object and
#' raises a \code{gdpar_information_error} warning class but does
#' \emph{not} abort the fit. The user is expected to inspect the
#' returned report and decide whether to discard the slot's anchor.
#'
#' @param fit A cmdstanr fit object (from \code{cs_model$sample()}).
#' @param family A \code{gdpar_family} object whose K-individual slots
#'   have been promoted to per-observation scope by
#'   \code{.gdpar_promote_scope_per_observation()}.
#' @param slot_names Character vector of length K with the canonical
#'   slot names.
#' @param use_groups Integer scalar (0 or 1) declaring whether the fit
#'   used per-group anchors (drives the choice between
#'   \code{theta_ref_k} and \code{mu_theta_ref_k}).
#' @param prior A \code{gdpar_prior} object; used to consult
#'   \code{priors_by_kind} overrides (currently inert; the helper falls
#'   back to canonical kinds).
#' @return List with \code{passed} (logical), \code{thresholds} (named
#'   numeric vector with \code{warn} = 0.5 and
#'   \code{information_error} = 0.1), \code{per_slot} (list of length
#'   K; each entry has \code{slot}, \code{var_post}, \code{var_prior},
#'   \code{contraction}, \code{status} in
#'   \code{c("pass", "warn", "information_error", "skipped")},
#'   \code{message}).
#' @keywords internal
#' @noRd
.compute_info_ratio_K <- function(fit, family, slot_names,
                                   use_groups, prior) {
  K <- length(slot_names)
  draws <- tryCatch(
    fit$draws(format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws)) {
    return(list(
      passed     = TRUE,
      thresholds = c(warn = 0.5, information_error = 0.1),
      per_slot   = lapply(seq_len(K), function(k) {
        list(slot = slot_names[k], var_post = NA_real_,
             var_prior = NA_real_, contraction = NA_real_,
             status = "skipped",
             message = "fit$draws() unavailable")
      })
    ))
  }
  param_root <- if (use_groups == 1L) "mu_theta_ref_k" else "theta_ref_k"
  param_var_names <- if (use_groups == 1L) {
    paste0(param_root, "[1,", seq_len(K), "]")
  } else {
    paste0(param_root, "[1,", seq_len(K), "]")
  }
  individual_specs <- family$param_specs[
    vapply(family$param_specs, function(s) {
      s$scope %in% c("per_observation", "per_group")
    }, logical(1L))
  ]
  if (length(individual_specs) != K) {
    individual_specs <- family$param_specs[seq_len(K)]
  }
  per_slot <- vector("list", K)
  any_warn <- FALSE
  any_info_err <- FALSE
  for (k in seq_len(K)) {
    pv <- param_var_names[k]
    spec_k <- individual_specs[[k]]
    kind <- spec_k$prior_canonical_kind
    var_prior <- .gdpar_canonical_prior_variance(kind)
    if (!is.finite(var_prior)) {
      per_slot[[k]] <- list(
        slot        = slot_names[k],
        var_post    = NA_real_,
        var_prior   = NA_real_,
        contraction = NA_real_,
        status      = "skipped",
        message     = sprintf(
          "no canonical prior variance for kind '%s'", kind
        )
      )
      next
    }
    draws_k <- if (pv %in% colnames(draws)) draws[, pv] else NULL
    if (is.null(draws_k) || length(draws_k) < 4L) {
      per_slot[[k]] <- list(
        slot        = slot_names[k],
        var_post    = NA_real_,
        var_prior   = var_prior,
        contraction = NA_real_,
        status      = "skipped",
        message     = sprintf("draws for '%s' unavailable", pv)
      )
      next
    }
    var_post <- stats::var(as.numeric(draws_k))
    contraction <- 1 - var_post / var_prior
    if (!is.finite(contraction)) {
      per_slot[[k]] <- list(
        slot        = slot_names[k],
        var_post    = var_post,
        var_prior   = var_prior,
        contraction = NA_real_,
        status      = "skipped",
        message     = "non-finite contraction"
      )
      next
    }
    status <- if (contraction < 0.1) {
      any_info_err <- TRUE
      "information_error"
    } else if (contraction < 0.5) {
      any_warn <- TRUE
      "warn"
    } else {
      "pass"
    }
    per_slot[[k]] <- list(
      slot        = slot_names[k],
      var_post    = var_post,
      var_prior   = var_prior,
      contraction = contraction,
      status      = status,
      message     = sprintf(
        "contraction = %.3f (kind '%s', var_prior = %.3f)",
        contraction, kind, var_prior
      )
    )
  }
  names(per_slot) <- slot_names
  list(
    passed     = !any_info_err && !any_warn,
    any_warn   = any_warn,
    any_info_error = any_info_err,
    thresholds = c(warn = 0.5, information_error = 0.1),
    per_slot   = per_slot
  )
}

#' Print method for gdpar_identifiability_report objects
#'
#' @param x An object of class \code{gdpar_identifiability_report}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_identifiability_report <- function(x, ...) {
  cat("<gdpar_identifiability_report>\n")
  cat("  passed             : ", x$passed, "\n", sep = "")
  if (!is.na(x$lambda_min)) {
    cat("  lambda_min         : ", format(x$lambda_min, digits = 4),
        "\n", sep = "")
    cat("  lambda_max         : ", format(x$lambda_max, digits = 4),
        "\n", sep = "")
    cat("  condition_number   : ", format(x$condition_number, digits = 4),
        "\n", sep = "")
    cat("  tol                : ", format(x$tol_used, digits = 4),
        "\n", sep = "")
  }
  if (!is.null(x$rigor_used)) {
    cat("  rigor              : ", x$rigor_used, "\n", sep = "")
  }
  cat("  ", x$message, "\n", sep = "")
  if (!x$passed && !is.null(x$collinear_directions)) {
    cat("  collinear directions (C1-C4):\n")
    for (d in x$collinear_directions) {
      cat("    eigenvalue ", format(d$eigenvalue, digits = 4),
          ":\n", sep = "")
      for (k in seq_along(d$columns)) {
        cat("      ", format(d$coefficients[k], digits = 3,
                              width = 7),
            " * ", d$columns[k], "\n", sep = "")
      }
    }
  }
  if (!is.null(x$c4_bis)) {
    cat("  C4-bis (per-coordinate cross-component):\n")
    for (pk in x$c4_bis$per_k) {
      cat("    coord ", pk$coord, " (rigor=", pk$rigor,
          ") passed=", pk$passed, sep = "")
      if (is.finite(pk$condition_number)) {
        cat(" cond=", format(pk$condition_number, digits = 4),
            sep = "")
      }
      if (length(pk$shared_cols) > 0L) {
        cat(" shared_cols={",
            paste(pk$shared_cols, collapse = ","), "}",
            sep = "")
      }
      cat("\n")
      if (!isTRUE(pk$passed) && !is.null(pk$collinear_directions)) {
        for (d in pk$collinear_directions) {
          cat("      eigenvalue ", format(d$eigenvalue, digits = 4),
              ":\n", sep = "")
          for (kk in seq_along(d$columns)) {
            cat("        ",
                format(d$coefficients[kk], digits = 3, width = 7),
                " * ", d$columns[kk], "\n", sep = "")
          }
        }
      }
    }
  }
  if (!is.null(x$did_pre_fit)) {
    cat("  D-ID pre-fit (Block 8 decision 3C layer ii):\n")
    cat("    K (individual params)  : ", x$did_pre_fit$K, "\n", sep = "")
    cat("    passed                 : ", x$did_pre_fit$passed, "\n",
        sep = "")
    for (pp in x$did_pre_fit$per_param) {
      cat("    param '", pp$name, "' (scope=", pp$scope,
          ", kind=", pp$prior_canonical_kind, "): ", pp$did_status,
          sep = "")
      if (!is.na(pp$did_condition)) {
        cat(" -- ", pp$did_condition, sep = "")
      }
      cat("\n")
    }
    if (!is.null(x$did_pre_fit$symbolic_separability)) {
      cat("    symbolic separability  : ",
          x$did_pre_fit$symbolic_separability$passed, "\n", sep = "")
      if (!isTRUE(x$did_pre_fit$symbolic_separability$passed)) {
        cat("      overlapping kinds: ",
            paste(x$did_pre_fit$symbolic_separability$overlapping_kinds,
                  collapse = ", "),
            "\n", sep = "")
      }
    }
  }
  invisible(x)
}

#' Internal: check (C7) group-anchor aliasing
#'
#' Implements the condition (C7) of Block 6.5: when the user activates
#' the per-group hierarchical anchor (group argument of gdpar() non-NULL),
#' neither Z_a nor Z_b may contain columns that are linearly dependent
#' with the group indicator. The check has two layers:
#'
#' 1. Constant-per-group columns: any column of Z_a or Z_b that is
#'    constant within every level of group_id is rank-deficient with
#'    the group dummy matrix and creates a perfect alias with
#'    theta_ref[g]. This catches the typical mistake of including
#'    factor(group) (or any deterministic function of the group)
#'    inside a or b.
#' 2. Joint rank check: rank(cbind(G, Z)) < ncol(G) + ncol(Z) (with G
#'    the model.matrix of factor(group_id) without intercept and Z the
#'    centered design block). Catches indirect aliases where no single
#'    column is constant per group but a non-trivial linear combination
#'    is.
#'
#' Both layers operate on the centered design returned by
#' build_amm_design (Z_a, Z_b already have their column means removed),
#' which is what the Stan template consumes.
#'
#' @param design List returned by build_amm_design (univariate) or
#'   .build_amm_design_multi (multivariate, with Z_a_list and
#'   Z_b_list slots).
#' @param group_id Integer vector with values in 1..J_groups.
#' @param group_var_name Character scalar with the user-facing name of
#'   the group variable, used for error messages.
#' @param tol Numeric tolerance for the rank check (eigenvalue ratio).
#' @return Invisible NULL when the check passes. Aborts with
#'   gdpar_input_error class on violation.
#' @keywords internal
#' @noRd
.check_group_aliasing_c7 <- function(design, group_id, group_var_name,
                                     tol = 1e-8) {
  if (is.null(group_id)) return(invisible(NULL))
  group_fac <- as.factor(group_id)
  J_groups <- nlevels(group_fac)
  if (J_groups < 2L) return(invisible(NULL))
  group_int <- as.integer(group_fac)
  has_multi_design <- !is.null(design$Z_a_list)
  if (has_multi_design) {
    Z_a_blocks <- design$Z_a_list
    Z_b_blocks <- design$Z_b_list
    Z_a_names_list <- design$Z_a_names_list
    Z_b_names_list <- design$Z_b_names_list
    p <- length(Z_a_blocks)
    for (k in seq_len(p)) {
      .check_c7_one_block(Z_a_blocks[[k]],
                          Z_a_names_list[[k]],
                          component = "a", coord = k,
                          group_int = group_int, J_groups = J_groups,
                          group_var_name = group_var_name, tol = tol)
      .check_c7_one_block(Z_b_blocks[[k]],
                          Z_b_names_list[[k]],
                          component = "b", coord = k,
                          group_int = group_int, J_groups = J_groups,
                          group_var_name = group_var_name, tol = tol)
    }
  } else {
    .check_c7_one_block(design$Z_a, design$Z_a_names,
                        component = "a", coord = NA_integer_,
                        group_int = group_int, J_groups = J_groups,
                        group_var_name = group_var_name, tol = tol)
    .check_c7_one_block(design$Z_b, design$Z_b_names,
                        component = "b", coord = NA_integer_,
                        group_int = group_int, J_groups = J_groups,
                        group_var_name = group_var_name, tol = tol)
  }
  invisible(NULL)
}

#' Internal: (C7) per-block check helper
#'
#' Applies the two-layer aliasing check to one design block (Z_a or Z_b,
#' either univariate or per-coord-k under the multi path).
#'
#' @keywords internal
#' @noRd
.check_c7_one_block <- function(Z, Z_names, component, coord,
                                group_int, J_groups, group_var_name,
                                tol) {
  if (is.null(Z) || ncol(Z) == 0L) return(invisible(NULL))
  n <- nrow(Z)
  if (n == 0L) return(invisible(NULL))
  within_var <- vapply(seq_len(ncol(Z)), function(j) {
    col <- Z[, j]
    sums <- tapply(col, group_int, function(v) {
      if (length(v) < 2L) 0 else stats::var(v)
    })
    max(sums, na.rm = TRUE)
  }, numeric(1L))
  constant_cols <- which(within_var <= tol)
  if (length(constant_cols) > 0L) {
    bad_names <- Z_names[constant_cols]
    coord_tag <- if (is.na(coord)) "" else sprintf(" (coord k = %d)",
                                                    coord)
    gdpar_abort(
      sprintf(
        paste(
          "Identifiability condition (C7) of Block 6.5 violated: the",
          "following column(s) of the %s-component design%s are",
          "constant within every level of group '%s' and therefore",
          "perfectly alias the per-group anchor theta_ref[g]. Remove",
          "them from the %s spec, or fit without the group argument: %s."
        ),
        component, coord_tag, group_var_name, component,
        paste(bad_names, collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(
        component = component,
        coord = coord,
        group_var = group_var_name,
        aliased_columns = bad_names
      )
    )
  }
  G <- stats::model.matrix(~ as.factor(group_int) + 0)
  M <- cbind(G, Z)
  col_norms <- sqrt(colSums(M * M))
  zero_idx <- which(col_norms == 0)
  if (length(zero_idx) > 0L) {
    col_norms[zero_idx] <- 1
  }
  M_norm <- sweep(M, 2L, col_norms, "/")
  qr_M <- qr(M_norm, tol = tol)
  if (qr_M$rank < ncol(M_norm)) {
    coord_tag <- if (is.na(coord)) "" else sprintf(" (coord k = %d)",
                                                    coord)
    gdpar_abort(
      sprintf(
        paste(
          "Identifiability condition (C7) of Block 6.5 violated: the",
          "joint matrix [G | Z_%s]%s is rank-deficient (rank %d of %d).",
          "This indicates an indirect alias between the per-group anchor",
          "and a linear combination of the %s columns. Remove the",
          "aliased columns from the %s spec, or fit without the group",
          "argument."
        ),
        component, coord_tag, qr_M$rank, ncol(M_norm), component,
        component
      ),
      class = "gdpar_input_error",
      data = list(
        component = component,
        coord = coord,
        group_var = group_var_name,
        rank = qr_M$rank,
        ncol = ncol(M_norm)
      )
    )
  }
  invisible(NULL)
}


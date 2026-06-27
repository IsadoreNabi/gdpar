#' Adapter for external meta-learners used by gdpar_compare_meta_learners
#'
#' Constructor of the pluggable contract through which the comparator
#' \code{\link{gdpar_compare_meta_learners}} dispatches to external
#' meta-learner implementations such as \pkg{grf} or EconML. The
#' constructor returns an object of class
#' \code{gdpar_meta_learner_adapter} that the comparator orchestrates;
#' two reference adapters distributed with the package
#' (\code{\link{gdpar_adapter_grf}} and \code{\link{gdpar_adapter_econml}})
#' are produced by calling this constructor with the appropriate
#' closures. Users can build additional adapters by passing their own
#' fitting and (optionally) prediction closures, which makes the
#' comparator open to integration with other meta-learner ecosystems
#' (e.g. \pkg{DoubleML}, custom doubly-robust estimators) without
#' modification of the package code.
#'
#' The contract is intentionally two-layered. The mandatory
#' \code{fit_predict_fun} carries the full fit-and-predict cycle and is
#' the only function required to build a valid adapter. The optional
#' \code{predict_fun} re-evaluates the meta-learner on a new evaluation
#' grid by reusing the fitted state returned by \code{fit_predict_fun};
#' when an adapter exposes \code{predict_fun}, the comparator method
#' \code{\link{predict.gdpar_meta_learner_comparison}} dispatches to it
#' and avoids a costly re-fit. When \code{predict_fun} is \code{NULL},
#' the comparator falls back to \code{fit_predict_fun} and emits a
#' \code{gdpar_diagnostic_warning} announcing the re-fit.
#'
#' @param name Character scalar with the adapter name. Must be
#'   non-empty and unique within a single call to
#'   \code{gdpar_compare_meta_learners}.
#' @param fit_predict_fun Function with signature
#'   \code{function(X, Y, T, X_newdata, level, seed_run)} that fits
#'   the meta-learner on \code{(X, Y, T)} and predicts the CATE on
#'   \code{X_newdata}. \code{X} and \code{X_newdata} are data frames
#'   of covariates (the adapter is responsible for any internal
#'   conversion to matrices or to language-native arrays such as
#'   \code{numpy.ndarray}); \code{Y} is a numeric vector and \code{T}
#'   a 0/1 integer vector. \code{level} is the nominal credible level
#'   inherited from the bridge (default 0.95); \code{seed_run} is the
#'   per-method seed propagated by the comparator. The function must
#'   return a list with components \code{cate_mean} (numeric vector),
#'   \code{cate_ci} (numeric matrix \code{[n_newdata, 2]} with columns
#'   \code{lower}, \code{upper}, or \code{NULL} when the adapter does
#'   not expose a native CI), \code{state} (opaque object cached for
#'   the optional \code{predict_fun}; may be \code{NULL}), and
#'   \code{notes} (character vector with method-specific diagnostics
#'   or warnings emitted during the fit).
#' @param predict_fun Optional function with signature
#'   \code{function(state, X_newdata, level)} that re-evaluates the
#'   meta-learner on \code{X_newdata} using the cached \code{state}
#'   produced by \code{fit_predict_fun}. Must return a list with
#'   components \code{cate_mean} and \code{cate_ci} (same shape as in
#'   \code{fit_predict_fun}). \code{NULL} (default) signals that the
#'   adapter does not support re-prediction; the comparator falls
#'   back to \code{fit_predict_fun} with a diagnostic warning.
#' @param requires_r Character vector of R packages that the adapter
#'   needs (checked via \code{requireNamespace} before the fit).
#'   Default \code{character(0)}.
#' @param requires_py Character vector of Python modules that the
#'   adapter needs (checked via
#'   \code{reticulate::py_module_available} before the fit, which
#'   itself requires the R package \pkg{reticulate} to be present).
#'   Default \code{character(0)}.
#' @param native_ci Logical scalar. \code{TRUE} when the adapter
#'   returns a native CI in \code{cate_ci}, \code{FALSE} when the
#'   adapter does not produce one (the comparator does not synthesize
#'   intervals in that case; the slot is left at \code{NULL}).
#' @param description Optional character scalar with a one-line
#'   human-readable description of the adapter, used by the
#'   \code{print} method.
#'
#' @return An object of class \code{gdpar_meta_learner_adapter} with
#'   components \code{name}, \code{fit_predict_fun},
#'   \code{predict_fun}, \code{requires_r}, \code{requires_py},
#'   \code{native_ci}, \code{description}.
#'
#' @seealso \code{\link{gdpar_compare_meta_learners}},
#'   \code{\link{gdpar_adapter_grf}}, \code{\link{gdpar_adapter_econml}}.
#'
#' @examples
#' fit_pred <- function(X, Y, T, X_newdata, level, seed_run) {
#'   m_t <- stats::lm(Y[T == 1L] ~ ., data = X[T == 1L, , drop = FALSE])
#'   m_c <- stats::lm(Y[T == 0L] ~ ., data = X[T == 0L, , drop = FALSE])
#'   p_t <- stats::predict(m_t, newdata = X_newdata)
#'   p_c <- stats::predict(m_c, newdata = X_newdata)
#'   list(cate_mean = as.numeric(p_t - p_c), cate_ci = NULL,
#'        state = list(m_t = m_t, m_c = m_c), notes = character(0))
#' }
#' lin_adapter <- gdpar_meta_learner_adapter(
#'   name = "lm_t_learner", fit_predict_fun = fit_pred,
#'   native_ci = FALSE,
#'   description = "Linear T-learner via stats::lm (example only)"
#' )
#' print(lin_adapter)
#'
#' @export
gdpar_meta_learner_adapter <- function(name,
                                       fit_predict_fun,
                                       predict_fun = NULL,
                                       requires_r = character(0L),
                                       requires_py = character(0L),
                                       native_ci = FALSE,
                                       description = NULL) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    gdpar_abort(
      "Argument 'name' must be a non-empty character scalar.",
      class = "gdpar_input_error",
      data = list(received = name)
    )
  }
  if (!is.function(fit_predict_fun)) {
    gdpar_abort(
      "Argument 'fit_predict_fun' must be a function.",
      class = "gdpar_input_error"
    )
  }
  if (!is.null(predict_fun) && !is.function(predict_fun)) {
    gdpar_abort(
      "Argument 'predict_fun' must be a function or NULL.",
      class = "gdpar_input_error"
    )
  }
  if (!is.character(requires_r)) {
    gdpar_abort(
      "Argument 'requires_r' must be a character vector.",
      class = "gdpar_input_error"
    )
  }
  if (!is.character(requires_py)) {
    gdpar_abort(
      "Argument 'requires_py' must be a character vector.",
      class = "gdpar_input_error"
    )
  }
  if (!is.logical(native_ci) || length(native_ci) != 1L ||
      is.na(native_ci)) {
    gdpar_abort(
      "Argument 'native_ci' must be a non-NA logical scalar.",
      class = "gdpar_input_error"
    )
  }
  if (!is.null(description) &&
      (!is.character(description) || length(description) != 1L)) {
    gdpar_abort(
      "Argument 'description' must be NULL or a character scalar.",
      class = "gdpar_input_error"
    )
  }
  obj <- list(
    name = name,
    fit_predict_fun = fit_predict_fun,
    predict_fun = predict_fun,
    requires_r = requires_r,
    requires_py = requires_py,
    native_ci = as.logical(native_ci),
    description = description
  )
  class(obj) <- c("gdpar_meta_learner_adapter", "list")
  obj
}

#' Test whether an object is a gdpar_meta_learner_adapter
#'
#' @param x Object to test.
#' @return \code{TRUE} when \code{x} inherits from class
#'   \code{gdpar_meta_learner_adapter}, \code{FALSE} otherwise.
#' @seealso \code{\link{gdpar_meta_learner_adapter}}.
#' @examples
#' a <- gdpar_meta_learner_adapter(
#'   name = "dummy",
#'   fit_predict_fun = function(X, Y, T, X_newdata, level, seed_run) {
#'     list(cate_mean = rep(0, nrow(X_newdata)), cate_ci = NULL,
#'          state = NULL, notes = character(0))
#'   }
#' )
#' is_gdpar_meta_learner_adapter(a)
#' is_gdpar_meta_learner_adapter(list())
#' @export
is_gdpar_meta_learner_adapter <- function(x) {
  inherits(x, "gdpar_meta_learner_adapter")
}

#' Print method for gdpar_meta_learner_adapter
#'
#' @param x A \code{gdpar_meta_learner_adapter} object.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_meta_learner_adapter <- function(x, ...) {
  cat("<gdpar_meta_learner_adapter>\n")
  cat("  name                : ", x$name, "\n", sep = "")
  if (!is.null(x$description)) {
    cat("  description         : ", x$description, "\n", sep = "")
  }
  cat("  requires (R)        : ",
      if (length(x$requires_r) == 0L) "<none>" else
        paste(x$requires_r, collapse = ", "),
      "\n", sep = "")
  cat("  requires (Python)   : ",
      if (length(x$requires_py) == 0L) "<none>" else
        paste(x$requires_py, collapse = ", "),
      "\n", sep = "")
  cat("  native CI           : ", x$native_ci, "\n", sep = "")
  cat("  supports predict    : ",
      !is.null(x$predict_fun), "\n", sep = "")
  invisible(x)
}

#' Internal: verify the adapter's R and Python requirements
#'
#' Returns a list with components \code{ok} (logical scalar),
#' \code{missing_r} (character vector of missing R packages), and
#' \code{missing_py} (character vector of missing Python modules).
#' Python availability is checked via
#' \code{reticulate::py_module_available}; if \pkg{reticulate} itself
#' is missing, all entries of \code{requires_py} are reported as
#' missing without attempting to load Python.
#'
#' @keywords internal
#' @noRd
.check_adapter_requirements <- function(adapter) {
  missing_r <- adapter$requires_r[
    !vapply(adapter$requires_r,
            function(pkg) requireNamespace(pkg, quietly = TRUE),
            logical(1L))
  ]
  if (length(adapter$requires_py) == 0L) {
    missing_py <- character(0L)
  } else if (!requireNamespace("reticulate", quietly = TRUE)) {
    missing_py <- adapter$requires_py
  } else {
    if (exists("py_require", envir = asNamespace("reticulate"),
               inherits = FALSE)) {
      tryCatch(reticulate::py_require(adapter$requires_py),
               error = function(e) NULL)
    }
    missing_py <- adapter$requires_py[
      !vapply(adapter$requires_py,
              function(mod) {
                tryCatch(
                  reticulate::py_module_available(mod),
                  error = function(e) FALSE
                )
              },
              logical(1L))
    ]
  }
  list(
    ok = (length(missing_r) == 0L) && (length(missing_py) == 0L),
    missing_r  = missing_r,
    missing_py = missing_py
  )
}

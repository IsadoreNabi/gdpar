#' Reference adapter: EconML CausalForestDML for gdpar_compare_meta_learners
#'
#' Build a \code{gdpar_meta_learner_adapter} that wraps a Python-side
#' EconML estimator (Chernozhukov et al., 2018) via \pkg{reticulate}
#' for use with \code{\link{gdpar_compare_meta_learners}}. The default
#' estimator is \code{econml.dml.CausalForestDML}, the orthogonal
#' double-machine-learning causal forest of Athey, Tibshirani, and
#' Wager (2019). The adapter exposes both the mandatory
#' \code{fit_predict_fun} and the optional \code{predict_fun}; the
#' latter reuses the fitted Python estimator on a fresh evaluation
#' grid without a refit. Native CIs are produced by EconML's
#' \code{effect_interval(X, alpha = 1 - level)} method.
#'
#' The Python module \code{econml} is in \code{Suggests} and must be
#' installed by the user outside of the package (e.g.
#' \code{reticulate::py_install("econml")} or manual installation in
#' the active Python environment). The adapter aborts cleanly when
#' \pkg{reticulate} or the \code{econml} module is unavailable.
#'
#' Cached state caveat: the returned \code{state} carries a reference
#' to a Python object. The reference is valid for the duration of the
#' R session in which the bridge was built; serializing the
#' comparison via \code{saveRDS} and reloading in a fresh R session
#' invalidates the Python reference and the \code{predict_fun} aborts
#' with \code{gdpar_unsupported_feature_error} when invoked on a
#' restored state. Rebuild the comparison in such cases.
#'
#' @param estimator Character scalar identifying the EconML estimator.
#'   Default \code{"CausalForestDML"}. Other identifiers can be added
#'   in future sub-phases.
#' @param n_estimators Integer scalar; number of trees in the EconML
#'   causal forest. Default \code{1000L}.
#' @param model_y Optional Python model object for the outcome stage
#'   (\code{econml}'s \code{model_y} argument). Default \code{NULL}
#'   delegates to EconML's default (a regression-tree-based model).
#' @param model_t Optional Python model object for the treatment
#'   stage (\code{econml}'s \code{model_t} argument). Default
#'   \code{NULL} delegates to EconML's default.
#' @param seed Optional integer scalar with the seed propagated to
#'   the EconML estimator (\code{random_state}). Default \code{NULL}.
#'
#' @return A \code{gdpar_meta_learner_adapter} object with
#'   \code{requires_r = "reticulate"},
#'   \code{requires_py = "econml"}, \code{native_ci = TRUE}, and both
#'   \code{fit_predict_fun} and \code{predict_fun} populated.
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen,
#' C., Newey, W., and Robins, J. (2018). Double/debiased machine
#' learning for treatment and structural parameters. \emph{The
#' Econometrics Journal}, 21(1), C1-C68.
#'
#' Athey, S., Tibshirani, J., and Wager, S. (2019). Generalized random
#' forests. \emph{The Annals of Statistics}, 47(2), 1148-1178.
#'
#' @seealso \code{\link{gdpar_compare_meta_learners}},
#'   \code{\link{gdpar_meta_learner_adapter}},
#'   \code{\link{gdpar_adapter_grf}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("reticulate", quietly = TRUE)) {
#'   adapter <- gdpar_adapter_econml(n_estimators = 200L)
#'   print(adapter)
#' }
#' }
#'
#' @export
gdpar_adapter_econml <- function(estimator = "CausalForestDML",
                                  n_estimators = 1000L,
                                  model_y = NULL,
                                  model_t = NULL,
                                  seed = NULL) {
  if (!is.character(estimator) || length(estimator) != 1L ||
      !nzchar(estimator)) {
    gdpar_abort("Argument 'estimator' must be a non-empty character scalar.",
                class = "gdpar_input_error")
  }
  if (!identical(estimator, "CausalForestDML")) {
    gdpar_abort(
      sprintf(
        "Estimator '%s' is not supported by gdpar_adapter_econml in Sub-phase 8.5.B. Use 'CausalForestDML' or define a custom adapter via gdpar_meta_learner_adapter().",
        estimator
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(estimator = estimator)
    )
  }
  assert_count(n_estimators, "n_estimators")
  if (!is.null(seed)) {
    assert_numeric_scalar(seed, "seed", lower = 1,
                          upper = .Machine$integer.max)
  }
  hp <- list(
    estimator = estimator,
    n_estimators = as.integer(n_estimators),
    model_y = model_y,
    model_t = model_t,
    seed = if (is.null(seed)) NULL else as.integer(seed)
  )
  fit_predict_fun <- function(X, Y, T, X_newdata, level, seed_run) {
    require_suggested("reticulate", "fit gdpar_adapter_econml")
    if (exists("py_require", envir = asNamespace("reticulate"),
               inherits = FALSE)) {
      tryCatch(reticulate::py_require("econml"),
               error = function(e) NULL)
    }
    if (!isTRUE(tryCatch(
      reticulate::py_module_available("econml"),
      error = function(e) FALSE
    ))) {
      gdpar_abort(
        paste(
          "Python module 'econml' is not available in the active",
          "reticulate Python environment. Install it via",
          "reticulate::py_install('econml') or manually, and rerun."
        ),
        class = "gdpar_missing_dependency_error",
        data = list(python_module = "econml")
      )
    }
    X_mat <- .econml_to_matrix(X)
    X_new_mat <- .econml_to_matrix(X_newdata,
                                    template = attr(X_mat, "template"))
    np <- reticulate::import("numpy", convert = FALSE)
    econml_dml <- reticulate::import("econml.dml", convert = FALSE)
    eff_seed <- if (!is.null(seed_run)) as.integer(seed_run)
                else hp$seed
    args <- list(n_estimators = hp$n_estimators)
    if (!is.null(eff_seed)) args$random_state <- eff_seed
    if (!is.null(hp$model_y)) args$model_y <- hp$model_y
    if (!is.null(hp$model_t)) args$model_t <- hp$model_t
    est <- do.call(econml_dml$CausalForestDML, args)
    est$fit(
      Y = np$asarray(as.numeric(Y)),
      T = np$asarray(as.numeric(T)),
      X = np$asarray(X_mat)
    )
    eff <- est$effect(np$asarray(X_new_mat))
    cate_mean <- as.numeric(reticulate::py_to_r(eff))
    alpha <- 1 - level
    ci_py <- est$effect_interval(np$asarray(X_new_mat), alpha = alpha)
    ci_list <- reticulate::py_to_r(ci_py)
    cate_ci <- cbind(lower = as.numeric(ci_list[[1L]]),
                     upper = as.numeric(ci_list[[2L]]))
    list(
      cate_mean = cate_mean,
      cate_ci   = cate_ci,
      state     = list(
        estimator = est,
        template  = attr(X_mat, "template")
      ),
      notes     = character(0L)
    )
  }
  predict_fun <- function(state, X_newdata, level) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      gdpar_abort(
        "Package 'reticulate' is required to reuse a cached EconML state.",
        class = "gdpar_missing_dependency_error",
        data = list(package = "reticulate")
      )
    }
    if (exists("py_require", envir = asNamespace("reticulate"),
               inherits = FALSE)) {
      tryCatch(reticulate::py_require("econml"),
               error = function(e) NULL)
    }
    if (is.null(state) || is.null(state$estimator)) {
      gdpar_abort(
        "Cached state for the econml adapter is empty; refit before predicting.",
        class = "gdpar_internal_error"
      )
    }
    if (!inherits(state$estimator, "python.builtin.object")) {
      gdpar_abort(
        paste(
          "Cached EconML estimator has lost its Python reference",
          "(typically after saveRDS round-trip). Rebuild the",
          "comparison from a fresh fit."
        ),
        class = "gdpar_unsupported_feature_error"
      )
    }
    X_new_mat <- .econml_to_matrix(X_newdata, template = state$template)
    np <- reticulate::import("numpy", convert = FALSE)
    eff <- state$estimator$effect(np$asarray(X_new_mat))
    cate_mean <- as.numeric(reticulate::py_to_r(eff))
    alpha <- 1 - level
    ci_py <- state$estimator$effect_interval(np$asarray(X_new_mat),
                                              alpha = alpha)
    ci_list <- reticulate::py_to_r(ci_py)
    cate_ci <- cbind(lower = as.numeric(ci_list[[1L]]),
                     upper = as.numeric(ci_list[[2L]]))
    list(cate_mean = cate_mean, cate_ci = cate_ci)
  }
  gdpar_meta_learner_adapter(
    name = "econml",
    fit_predict_fun = fit_predict_fun,
    predict_fun = predict_fun,
    requires_r  = "reticulate",
    requires_py = "econml",
    native_ci   = TRUE,
    description = paste0(
      "econml.dml.", estimator, " (n_estimators = ", hp$n_estimators,
      ") via reticulate; native CIs from effect_interval()."
    )
  )
}

#' Internal: convert covariates to numeric matrix for EconML
#'
#' Mirror of \code{.grf_to_matrix} but kept separate so the two
#' adapters remain self-contained; both produce a fully numeric design
#' matrix and record a template so newdata aligns.
#'
#' @keywords internal
#' @noRd
.econml_to_matrix <- function(df, template = NULL) {
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    if (is.character(df[[j]])) df[[j]] <- as.factor(df[[j]])
  }
  if (is.null(template)) {
    if (ncol(df) == 0L) {
      gdpar_abort(
        "gdpar_adapter_econml requires at least one covariate; received a 0-column data frame.",
        class = "gdpar_input_error"
      )
    }
    mm <- stats::model.matrix(~ . - 1, data = df)
    template <- list(
      colnames = colnames(mm),
      factor_levels = lapply(df, function(col) {
        if (is.factor(col)) levels(col) else NULL
      })
    )
  } else {
    for (j in names(template$factor_levels)) {
      if (!is.null(template$factor_levels[[j]])) {
        if (!j %in% colnames(df)) {
          gdpar_abort(
            sprintf("Covariate '%s' missing from newdata for the econml adapter.", j),
            class = "gdpar_input_error",
            data = list(missing = j)
          )
        }
        df[[j]] <- factor(df[[j]], levels = template$factor_levels[[j]])
      }
    }
    mm <- stats::model.matrix(~ . - 1, data = df)
    if (!setequal(colnames(mm), template$colnames)) {
      missing_cols <- setdiff(template$colnames, colnames(mm))
      extra_cols   <- setdiff(colnames(mm), template$colnames)
      gdpar_abort(
        sprintf(
          "econml adapter: newdata produced an incompatible design (missing: %s; extra: %s).",
          if (length(missing_cols) == 0L) "<none>"
            else paste(missing_cols, collapse = ", "),
          if (length(extra_cols)   == 0L) "<none>"
            else paste(extra_cols,   collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(missing = missing_cols, extra = extra_cols)
      )
    }
    mm <- mm[, template$colnames, drop = FALSE]
  }
  attr(mm, "template") <- template
  mm
}

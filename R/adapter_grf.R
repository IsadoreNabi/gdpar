#' Reference adapter: grf::causal_forest for gdpar_compare_meta_learners
#'
#' Build a \code{gdpar_meta_learner_adapter} that wraps the R-side
#' causal forest of \pkg{grf} (Athey, Tibshirani, and Wager, 2019) for
#' use with \code{\link{gdpar_compare_meta_learners}}. The adapter
#' exposes both the mandatory \code{fit_predict_fun} and the optional
#' \code{predict_fun} so the comparator's \code{predict} method can
#' reuse the fitted causal forest on a fresh evaluation grid without
#' a refit. Native CIs are obtained by the normal approximation
#' \eqn{\widehat\tau(x) \pm z_{1-\alpha/2} \cdot \sqrt{\widehat{\mathrm{Var}}(\widehat\tau(x))}}
#' using \code{grf}'s built-in variance estimator
#' (\code{predict(..., estimate.variance = TRUE)}).
#'
#' Categorical covariates are not handled by \pkg{grf} directly; the
#' adapter coerces character columns to factors and then applies
#' \code{stats::model.matrix(~ . - 1, ...)} to obtain a fully numeric
#' design matrix. Numeric or factor inputs pass through unchanged.
#'
#' @param num_trees Integer scalar; number of trees in the forest.
#'   Default \code{2000L}, matching \pkg{grf}'s default.
#' @param sample_fraction Numeric scalar in \code{(0, 0.5]}; fraction
#'   of the training sample drawn for each tree. Default \code{0.5}.
#' @param mtry Optional integer scalar with the number of candidate
#'   variables per split; default \code{NULL} delegates to \pkg{grf}'s
#'   own default (\code{min(ceiling(sqrt(p) + 20), p)}).
#' @param honesty Logical scalar; whether to use honest splitting.
#'   Default \code{TRUE} (recommended; the \pkg{grf} confidence
#'   intervals are valid only under honesty).
#' @param seed Optional integer scalar with the seed propagated to
#'   \pkg{grf}'s internal RNG when the comparator's \code{seed_run} is
#'   \code{NULL}. Default \code{NULL}.
#'
#' @return A \code{gdpar_meta_learner_adapter} object with
#'   \code{requires_r = "grf"}, \code{native_ci = TRUE}, and both
#'   \code{fit_predict_fun} and \code{predict_fun} populated.
#'
#' @references
#' Athey, S., Tibshirani, J., and Wager, S. (2019). Generalized random
#' forests. \emph{The Annals of Statistics}, 47(2), 1148-1178.
#'
#' Wager, S., and Athey, S. (2018). Estimation and inference of
#' heterogeneous treatment effects using random forests.
#' \emph{Journal of the American Statistical Association}, 113(523),
#' 1228-1242.
#'
#' @seealso \code{\link{gdpar_compare_meta_learners}},
#'   \code{\link{gdpar_meta_learner_adapter}},
#'   \code{\link{gdpar_adapter_econml}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("grf", quietly = TRUE)) {
#'   adapter <- gdpar_adapter_grf(num_trees = 500L)
#'   print(adapter)
#' }
#' }
#'
#' @export
gdpar_adapter_grf <- function(num_trees = 2000L,
                              sample_fraction = 0.5,
                              mtry = NULL,
                              honesty = TRUE,
                              seed = NULL) {
  assert_count(num_trees, "num_trees")
  assert_numeric_scalar(sample_fraction, "sample_fraction",
                        lower = 1e-3, upper = 0.5)
  if (!is.null(mtry)) assert_count(mtry, "mtry")
  if (!is.logical(honesty) || length(honesty) != 1L || is.na(honesty)) {
    gdpar_abort("Argument 'honesty' must be a non-NA logical scalar.",
                class = "gdpar_input_error")
  }
  if (!is.null(seed)) {
    assert_numeric_scalar(seed, "seed", lower = 1,
                          upper = .Machine$integer.max)
  }
  hp <- list(num_trees = as.integer(num_trees),
             sample_fraction = as.numeric(sample_fraction),
             mtry = if (is.null(mtry)) NULL else as.integer(mtry),
             honesty = as.logical(honesty),
             seed = if (is.null(seed)) NULL else as.integer(seed))
  fit_predict_fun <- function(X, Y, T, X_newdata, level, seed_run) {
    require_suggested("grf", "fit gdpar_adapter_grf")
    X_mat <- .grf_to_matrix(X)
    X_new_mat <- .grf_to_matrix(X_newdata, template = attr(X_mat, "template"))
    eff_seed <- if (!is.null(seed_run)) as.integer(seed_run)
                else hp$seed
    args <- list(
      X = X_mat, Y = as.numeric(Y), W = as.numeric(T),
      num.trees = hp$num_trees,
      sample.fraction = hp$sample_fraction,
      honesty = hp$honesty
    )
    if (!is.null(hp$mtry)) args$mtry <- hp$mtry
    if (!is.null(eff_seed)) args$seed <- eff_seed
    cf <- do.call(grf::causal_forest, args)
    pred <- stats::predict(cf, newdata = X_new_mat,
                            estimate.variance = TRUE)
    cate_mean <- as.numeric(pred$predictions)
    z <- stats::qnorm(1 - (1 - level) / 2)
    se <- sqrt(pmax(as.numeric(pred$variance.estimates), 0))
    cate_ci <- cbind(lower = cate_mean - z * se,
                     upper = cate_mean + z * se)
    list(
      cate_mean = cate_mean,
      cate_ci   = cate_ci,
      state     = list(
        forest   = cf,
        template = attr(X_mat, "template")
      ),
      notes     = character(0L)
    )
  }
  predict_fun <- function(state, X_newdata, level) {
    if (!requireNamespace("grf", quietly = TRUE)) {
      gdpar_abort(
        "Package 'grf' is required to reuse a cached causal_forest state.",
        class = "gdpar_missing_dependency_error",
        data = list(package = "grf")
      )
    }
    if (is.null(state) || is.null(state$forest)) {
      gdpar_abort(
        "Cached state for the grf adapter is empty; refit before predicting.",
        class = "gdpar_internal_error"
      )
    }
    X_new_mat <- .grf_to_matrix(X_newdata, template = state$template)
    pred <- stats::predict(state$forest, newdata = X_new_mat,
                            estimate.variance = TRUE)
    cate_mean <- as.numeric(pred$predictions)
    z <- stats::qnorm(1 - (1 - level) / 2)
    se <- sqrt(pmax(as.numeric(pred$variance.estimates), 0))
    cate_ci <- cbind(lower = cate_mean - z * se,
                     upper = cate_mean + z * se)
    list(cate_mean = cate_mean, cate_ci = cate_ci)
  }
  gdpar_meta_learner_adapter(
    name = "grf",
    fit_predict_fun = fit_predict_fun,
    predict_fun = predict_fun,
    requires_r = "grf",
    native_ci = TRUE,
    description = paste0(
      "grf::causal_forest (num_trees = ", hp$num_trees,
      ", honesty = ", hp$honesty,
      ") with normal-approximation CIs from estimate.variance."
    )
  )
}

#' Internal: convert a covariate data frame into a numeric matrix
#'
#' Characters are coerced to factors, factors are expanded via
#' \code{model.matrix(~ . - 1)}, and numeric columns pass through. A
#' \code{template} attribute records the column structure of the first
#' call so subsequent calls on \code{X_newdata} align identically;
#' aborts when a new column appears or a level is missing.
#'
#' @keywords internal
#' @noRd
.grf_to_matrix <- function(df, template = NULL) {
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    if (is.character(df[[j]])) df[[j]] <- as.factor(df[[j]])
  }
  if (is.null(template)) {
    if (ncol(df) == 0L) {
      gdpar_abort(
        "gdpar_adapter_grf requires at least one covariate; received a 0-column data frame.",
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
            sprintf("Covariate '%s' missing from newdata for the grf adapter.", j),
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
          "grf adapter: newdata produced an incompatible design (missing: %s; extra: %s).",
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

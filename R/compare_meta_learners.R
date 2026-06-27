#' Compare the AMM-side T-learner against external meta-learners
#'
#' Evaluate a fitted \code{\link{gdpar_causal_bridge}} object against a
#' user-supplied set of external meta-learner adapters (e.g.
#' \code{\link{gdpar_adapter_grf}} for \pkg{grf} R-side,
#' \code{\link{gdpar_adapter_econml}} for EconML Python-side) on a
#' common evaluation grid. The function does not refit either of the
#' two gdpar fits embedded in \code{bridge}; it only consumes the
#' bridge's CATE estimates and reconstructs the (X, T, Y) dataset
#' needed by the external adapters from the captured calls of the two
#' fits, or from an explicit \code{data} argument when the captured
#' calls cannot be resolved.
#'
#' The comparator is descriptive: it reports per-method posterior /
#' point CATE estimates together with their native CIs (when the
#' adapter exposes one) and three concordance metrics (RMSE, Pearson
#' correlation, mean absolute discrepancy) between every ordered pair
#' of methods over \code{cate_mean}. Tests of hypothesis and claims of
#' algorithmic equivalence are deliberately out of scope (the
#' inferential origin of each method differs); the interpretation of
#' the discrepancy is left to the user.
#'
#' @param bridge An object of class \code{gdpar_causal_bridge}
#'   produced by \code{\link{gdpar_causal_bridge}}.
#' @param methods A non-empty named or unnamed list of objects of
#'   class \code{gdpar_meta_learner_adapter}; the comparator orchestrates
#'   each adapter in turn. When the list is unnamed, the names of the
#'   methods are taken from the \code{name} field of each adapter.
#' @param newdata Optional data frame on which the CATE is evaluated.
#'   When \code{NULL} (default), the function reuses
#'   \code{bridge$newdata}.
#' @param data Optional list with components \code{X} (data frame of
#'   covariates), \code{T} (integer 0/1 vector), \code{Y} (numeric
#'   vector). Used when the captured calls of the two gdpar fits
#'   cannot be evaluated in the caller's environment (e.g. when the
#'   comparator is invoked from a wrapper that loses the data scope).
#'   When \code{NULL} (default), the function attempts to recover the
#'   training data of each arm from the bridge's stored fits.
#' @param seed Optional integer scalar propagated to each adapter as
#'   \code{seed_run}; default \code{NULL} leaves the RNG state to each
#'   adapter's discretion.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return An object of class \code{gdpar_meta_learner_comparison}
#'   with components \code{bridge_cate}, \code{external},
#'   \code{comparison}, \code{newdata}, \code{level}, \code{n_obs},
#'   \code{n_methods}, \code{call}, \code{meta}. See
#'   \code{\link{print.gdpar_meta_learner_comparison}} and
#'   \code{\link{summary.gdpar_meta_learner_comparison}} for the
#'   companion S3 methods.
#'
#' @section Scalar-outcome restriction:
#' The current scope of Sub-phase 8.5.B supports scalar outcomes only.
#' Bridges constructed from fits with \code{K > 1} (distributional
#' regression) or \code{p > 1} (multivariate response) are rejected
#' with \code{gdpar_unsupported_feature_error}; multi-output external
#' adapters are queued for Block 9 (see vignette
#' \code{v08c_meta_learner_comparison}, section "Limits").
#'
#' @section Dataset reconstruction:
#' When \code{data} is \code{NULL}, the helper
#' \code{.assemble_bridge_dataset} recovers the training data of each
#' arm via \code{eval(fit$call$data, eval_env)} (same mechanism used
#' by the bridge constructor in \code{\link{gdpar_causal_bridge}}),
#' identifies the outcome via the LHS of \code{fit$call$formula}, and
#' assembles a single (X, T, Y) dataset with \code{T = 1L} for the
#' treatment arm and \code{T = 0L} for the control arm. When the
#' evaluations fail, the helper aborts with \code{gdpar_input_error}
#' and instructs the user to pass \code{data} explicitly.
#'
#' @section Concordance metrics:
#' For every ordered pair of methods \eqn{(i, j)} (including the bridge
#' as a method indexed by \code{"bridge"}), the comparator reports
#' \deqn{\mathrm{RMSE}_{ij} = \sqrt{\mathrm{mean}((\widehat\tau_i - \widehat\tau_j)^2)}}
#' \deqn{\mathrm{Pearson}_{ij} = \mathrm{cor}(\widehat\tau_i, \widehat\tau_j)}
#' \deqn{\mathrm{MAD}_{ij} = \mathrm{mean}(\mid \widehat\tau_i - \widehat\tau_j \mid)}
#' computed on \code{cate_mean} only; CIs are not pooled across
#' methods because the inferential origin of each CI is heterogeneous
#' (posterior vs. asymptotic vs. bootstrap; see Appendix B of the
#' bridge vignette).
#'
#' @references
#' Kuenzel, S. R., Sekhon, J. S., Bickel, P. J., and Yu, B. (2019).
#' Metalearners for estimating heterogeneous treatment effects using
#' machine learning. \emph{Proceedings of the National Academy of
#' Sciences}, 116(10), 4156-4165.
#'
#' Athey, S., and Wager, S. (2019). Estimating treatment effects with
#' causal forests: An application. \emph{Observational Studies}, 5,
#' 37-51.
#'
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen,
#' C., Newey, W., and Robins, J. (2018). Double/debiased machine
#' learning for treatment and structural parameters. \emph{The
#' Econometrics Journal}, 21(1), C1-C68.
#'
#' @seealso \code{\link{gdpar_causal_bridge}},
#'   \code{\link{gdpar_meta_learner_adapter}},
#'   \code{\link{gdpar_adapter_grf}},
#'   \code{\link{gdpar_adapter_econml}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     requireNamespace("grf", quietly = TRUE)) {
#'   n <- 300L
#'   df <- data.frame(x1 = rnorm(2L * n))
#'   df$arm <- rep(c("treat", "ctrl"), each = n)
#'   df$y <- with(df, ifelse(arm == "treat", 0.5, 0) + 0.8 * x1 +
#'                    rnorm(2L * n, sd = 0.5))
#'   df_t <- subset(df, arm == "treat", select = -arm)
#'   df_c <- subset(df, arm == "ctrl",  select = -arm)
#'   fit_t <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df_t,
#'                  iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   fit_c <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df_c,
#'                  iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   bridge <- gdpar_causal_bridge(fit_t, fit_c,
#'                                 newdata = data.frame(x1 = seq(-2, 2, 0.2)))
#'   cmp <- gdpar_compare_meta_learners(bridge,
#'                                      methods = list(gdpar_adapter_grf()))
#'   print(cmp); summary(cmp)
#' }
#' }
#'
#' @export
gdpar_compare_meta_learners <- function(bridge, methods,
                                        newdata = NULL,
                                        data = NULL,
                                        seed = NULL,
                                        ...) {
  assert_inherits(bridge, "gdpar_causal_bridge", "bridge")
  .guard_scalar_outcome(bridge)
  if (!is.list(methods) || length(methods) == 0L) {
    gdpar_abort(
      "Argument 'methods' must be a non-empty list of gdpar_meta_learner_adapter objects.",
      class = "gdpar_input_error"
    )
  }
  bad <- vapply(methods, function(m) !is_gdpar_meta_learner_adapter(m),
                logical(1L))
  if (any(bad)) {
    gdpar_abort(
      sprintf(
        "All elements of 'methods' must be gdpar_meta_learner_adapter; found %d invalid element(s) at position(s) %s.",
        sum(bad), paste(which(bad), collapse = ", ")
      ),
      class = "gdpar_input_error"
    )
  }
  method_names <- names(methods)
  if (is.null(method_names) || any(!nzchar(method_names))) {
    method_names <- vapply(methods, function(m) m$name, character(1L))
  }
  if (anyDuplicated(method_names) > 0L) {
    gdpar_abort(
      sprintf(
        "Method names must be unique; received: (%s).",
        paste(method_names, collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(names = method_names)
    )
  }
  if (!is.null(seed)) {
    assert_numeric_scalar(seed, "seed", lower = 1, upper = .Machine$integer.max)
  }
  if (is.null(newdata)) {
    newdata_resolved <- bridge$newdata
  } else {
    assert_data_frame(newdata, "newdata")
    newdata_resolved <- newdata
  }
  ds <- .assemble_bridge_dataset(bridge, newdata = newdata_resolved,
                                  data = data,
                                  eval_env = parent.frame())
  level <- bridge$level
  bridge_cate_mean <- as.numeric(bridge$cate_mean)
  if (length(bridge_cate_mean) != nrow(newdata_resolved)) {
    if (is.null(newdata)) {
      gdpar_abort(
        "Internal error: bridge$cate_mean and bridge$newdata have inconsistent length.",
        class = "gdpar_internal_error"
      )
    }
    bridge_pred <- stats::predict(bridge, newdata = newdata_resolved,
                                  level = level, summary = "mean_ci")
    bridge_cate_mean <- as.numeric(bridge_pred$cate_mean)
    bridge_cate_ci <- bridge_pred$cate_ci
  } else {
    bridge_cate_ci <- bridge$cate_ci
  }
  external <- vector("list", length(methods))
  names(external) <- method_names
  for (k in seq_along(methods)) {
    adapter <- methods[[k]]
    nm <- method_names[k]
    req <- .check_adapter_requirements(adapter)
    if (!req$ok) {
      msg_r <- if (length(req$missing_r) > 0L) {
        sprintf("missing R package(s): %s",
                paste(req$missing_r, collapse = ", "))
      } else NULL
      msg_py <- if (length(req$missing_py) > 0L) {
        sprintf("missing Python module(s): %s",
                paste(req$missing_py, collapse = ", "))
      } else NULL
      gdpar_abort(
        sprintf(
          "Adapter '%s' has unmet requirements (%s). Install them and rerun.",
          nm, paste(c(msg_r, msg_py), collapse = "; ")
        ),
        class = "gdpar_missing_dependency_error",
        data = list(method = nm,
                    missing_r = req$missing_r,
                    missing_py = req$missing_py)
      )
    }
    t0 <- proc.time()
    result <- adapter$fit_predict_fun(
      X = ds$X, Y = ds$Y, T = ds$T,
      X_newdata = ds$X_newdata,
      level = level, seed_run = seed
    )
    t_elapsed <- (proc.time() - t0)[["elapsed"]]
    .validate_adapter_output(result, n_newdata = nrow(ds$X_newdata),
                             adapter_name = nm)
    external[[k]] <- list(
      cate_mean  = as.numeric(result$cate_mean),
      cate_ci    = result$cate_ci,
      method     = nm,
      native_ci  = adapter$native_ci && !is.null(result$cate_ci),
      time_sec   = unname(t_elapsed),
      notes      = if (is.null(result$notes)) character(0L)
                   else as.character(result$notes),
      state           = result$state,
      predict_fun     = adapter$predict_fun,
      fit_predict_fun = adapter$fit_predict_fun,
      has_predict_fun = !is.null(adapter$predict_fun)
    )
  }
  cate_list <- c(
    list(bridge = bridge_cate_mean),
    lapply(external, function(e) e$cate_mean)
  )
  comparison <- .compute_comparison_metrics(cate_list)
  obj <- list(
    bridge_cate = list(cate_mean = bridge_cate_mean,
                       cate_ci   = bridge_cate_ci),
    bridge      = bridge,
    external    = external,
    comparison  = comparison,
    newdata     = newdata_resolved,
    level       = level,
    n_obs       = nrow(newdata_resolved),
    n_methods   = length(methods),
    call        = match.call(),
    meta        = list(
      package_version = utils::packageVersion("gdpar"),
      date            = Sys.time(),
      seed            = seed,
      bridge_call     = bridge$call,
      adapters        = lapply(methods, function(m)
        list(name = m$name, requires_r = m$requires_r,
             requires_py = m$requires_py, native_ci = m$native_ci,
             has_predict_fun = !is.null(m$predict_fun)))
    )
  )
  class(obj) <- c("gdpar_meta_learner_comparison", "list")
  obj
}

#' Internal: assert the bridge is built from scalar-outcome fits
#' @keywords internal
#' @noRd
.guard_scalar_outcome <- function(bridge) {
  dim_kind <- bridge$meta$dim_kind %||% "scalar"
  if (!identical(dim_kind, "scalar")) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_compare_meta_learners (Sub-phase 8.5.B) supports scalar",
          "outcomes only; received bridge with dim_kind = '%s'.",
          "Multi-output external adapters (K > 1 or p > 1) are queued",
          "for Block 9."
        ),
        dim_kind
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(dim_kind = dim_kind)
    )
  }
  invisible(NULL)
}

#' Internal: assemble (X, T, Y, X_newdata) for external meta-learners
#'
#' Recovers the per-arm training data via \code{eval(fit$call$data, eval_env)}
#' when \code{data} is \code{NULL}, identifies the outcome from
#' \code{fit$call$formula}, and builds a single dataset stacking the
#' two arms with a treatment indicator. Returns a list with components
#' \code{X} (data frame), \code{T} (integer 0/1), \code{Y} (numeric),
#' \code{X_newdata} (data frame), \code{outcome_name} (character).
#' Aborts with \code{gdpar_input_error} when recovery fails or when
#' the user-supplied \code{data} is malformed.
#'
#' @keywords internal
#' @noRd
.assemble_bridge_dataset <- function(bridge, newdata, data, eval_env) {
  fit_t <- bridge$fits$treat
  fit_c <- bridge$fits$ctrl
  outcome_name <- .bridge_outcome_name(fit_t, fit_c)
  if (!is.null(data)) {
    if (!is.list(data) ||
        !all(c("X", "T", "Y") %in% names(data))) {
      gdpar_abort(
        "Argument 'data' must be a named list with components 'X', 'T', 'Y'.",
        class = "gdpar_input_error"
      )
    }
    X <- data$X
    if (!is.data.frame(X)) X <- as.data.frame(X, stringsAsFactors = FALSE)
    T_vec <- as.integer(data$T)
    Y_vec <- as.numeric(data$Y)
    if (length(T_vec) != length(Y_vec) || nrow(X) != length(Y_vec)) {
      gdpar_abort(
        sprintf(
          "Argument 'data': X, T, Y must share length; nrow(X) = %d, length(T) = %d, length(Y) = %d.",
          nrow(X), length(T_vec), length(Y_vec)
        ),
        class = "gdpar_input_error"
      )
    }
    if (!all(T_vec %in% c(0L, 1L))) {
      gdpar_abort(
        "Argument 'data$T' must contain only 0L and 1L.",
        class = "gdpar_input_error"
      )
    }
    X_newdata <- if (!is.null(data$X_newdata)) {
      if (!is.data.frame(data$X_newdata)) {
        as.data.frame(data$X_newdata, stringsAsFactors = FALSE)
      } else {
        data$X_newdata
      }
    } else {
      .extract_covariates(newdata, outcome_name)
    }
    return(list(X = X, T = T_vec, Y = Y_vec,
                X_newdata = X_newdata,
                outcome_name = outcome_name))
  }
  recover <- function(fit) {
    cl <- fit$call
    if (is.null(cl) || is.null(cl$data)) return(NULL)
    tryCatch(eval(cl$data, envir = eval_env), error = function(e) NULL)
  }
  data_t <- recover(fit_t)
  data_c <- recover(fit_c)
  if (is.null(data_t) || is.null(data_c) ||
      !is.data.frame(data_t) || !is.data.frame(data_c)) {
    gdpar_abort(
      paste(
        "gdpar_compare_meta_learners cannot recover the training data",
        "of one or both fits from the captured calls. Pass an explicit",
        "'data = list(X = ..., T = ..., Y = ..., X_newdata = ...)'",
        "argument."
      ),
      class = "gdpar_input_error",
      data = list(treat_recovered = !is.null(data_t),
                  ctrl_recovered  = !is.null(data_c))
    )
  }
  if (!(outcome_name %in% colnames(data_t)) ||
      !(outcome_name %in% colnames(data_c))) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' is missing from the recovered data of one or both arms.",
        outcome_name
      ),
      class = "gdpar_input_error",
      data = list(outcome_name = outcome_name)
    )
  }
  if (!identical(sort(colnames(data_t)), sort(colnames(data_c)))) {
    gdpar_abort(
      sprintf(
        paste(
          "Recovered training data of the two arms have non-matching",
          "column structures (treat = (%s), ctrl = (%s)). Pass an",
          "explicit 'data = list(X = ..., T = ..., Y = ..., X_newdata = ...)'."
        ),
        paste(colnames(data_t), collapse = ", "),
        paste(colnames(data_c), collapse = ", ")
      ),
      class = "gdpar_input_error"
    )
  }
  common_cols <- colnames(data_t)
  data_t <- data_t[, common_cols, drop = FALSE]
  data_c <- data_c[, common_cols, drop = FALSE]
  Y_vec <- c(as.numeric(data_t[[outcome_name]]),
             as.numeric(data_c[[outcome_name]]))
  T_vec <- c(rep(1L, nrow(data_t)), rep(0L, nrow(data_c)))
  X_full <- rbind(data_t, data_c)
  X <- .extract_covariates(X_full, outcome_name)
  X_newdata <- .extract_covariates(newdata, outcome_name)
  missing_cov <- setdiff(colnames(X), colnames(X_newdata))
  if (length(missing_cov) > 0L) {
    gdpar_abort(
      sprintf(
        "Newdata is missing covariate column(s): %s.",
        paste(missing_cov, collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(missing = missing_cov)
    )
  }
  X_newdata <- X_newdata[, colnames(X), drop = FALSE]
  list(X = X, T = T_vec, Y = Y_vec,
       X_newdata = X_newdata,
       outcome_name = outcome_name)
}

#' Internal: extract outcome name from the fit calls of a bridge
#' @keywords internal
#' @noRd
.bridge_outcome_name <- function(fit_t, fit_c) {
  pick <- function(fit) {
    cl <- fit$call
    if (is.null(cl) || is.null(cl$formula)) return(NA_character_)
    fm <- tryCatch(eval(cl$formula, envir = environment(fit$call)),
                    error = function(e) NULL)
    if (inherits(fm, "formula") && length(fm) == 3L) {
      return(as.character(fm[[2L]]))
    }
    if (inherits(cl$formula, "call") || inherits(cl$formula, "name")) {
      fm2 <- tryCatch(eval(cl$formula), error = function(e) NULL)
      if (inherits(fm2, "formula") && length(fm2) == 3L) {
        return(as.character(fm2[[2L]]))
      }
    }
    if (is.call(cl$formula) && length(cl$formula) == 3L &&
        identical(cl$formula[[1L]], as.name("~"))) {
      return(as.character(cl$formula[[2L]]))
    }
    NA_character_
  }
  n_t <- pick(fit_t)
  n_c <- pick(fit_c)
  if (is.na(n_t) && is.na(n_c)) {
    gdpar_abort(
      paste(
        "Cannot infer the outcome variable from the two gdpar fits;",
        "the captured formulas are not resolvable. Pass an explicit",
        "'data = list(X = ..., T = ..., Y = ..., X_newdata = ...)'."
      ),
      class = "gdpar_input_error"
    )
  }
  if (!is.na(n_t) && !is.na(n_c) && !identical(n_t, n_c)) {
    gdpar_abort(
      sprintf(
        "The two gdpar fits have non-matching outcome names ('%s' vs '%s').",
        n_t, n_c
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(outcome_treat = n_t, outcome_ctrl = n_c)
    )
  }
  if (is.na(n_t)) n_c else n_t
}

#' Internal: drop the outcome column and return covariates as data frame
#' @keywords internal
#' @noRd
.extract_covariates <- function(df, outcome_name) {
  if (!is.data.frame(df)) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  keep <- setdiff(colnames(df), outcome_name)
  df[, keep, drop = FALSE]
}

#' Internal: validate the shape of an adapter's output
#' @keywords internal
#' @noRd
.validate_adapter_output <- function(result, n_newdata, adapter_name) {
  if (!is.list(result) ||
      !("cate_mean" %in% names(result))) {
    gdpar_abort(
      sprintf(
        "Adapter '%s' did not return a list with a 'cate_mean' component.",
        adapter_name
      ),
      class = "gdpar_internal_error"
    )
  }
  if (length(result$cate_mean) != n_newdata ||
      !is.numeric(result$cate_mean)) {
    gdpar_abort(
      sprintf(
        "Adapter '%s' returned cate_mean of length %d; expected %d numeric values.",
        adapter_name, length(result$cate_mean), n_newdata
      ),
      class = "gdpar_internal_error"
    )
  }
  if (!is.null(result$cate_ci)) {
    if (!is.matrix(result$cate_ci) ||
        nrow(result$cate_ci) != n_newdata ||
        ncol(result$cate_ci) != 2L) {
      gdpar_abort(
        sprintf(
          "Adapter '%s' returned cate_ci with shape %s; expected %d-by-2 matrix.",
          adapter_name,
          paste(dim(result$cate_ci), collapse = " x "),
          n_newdata
        ),
        class = "gdpar_internal_error"
      )
    }
  }
  invisible(NULL)
}

#' Internal: compute the three concordance matrices over cate_list
#'
#' Returns a list with components \code{rmse}, \code{pearson},
#' \code{mad}, each an \code{m}-by-\code{m} matrix of pairwise metrics
#' over the \code{m = length(cate_list)} cate vectors. Diagonal entries
#' are 0 (or 1 for Pearson). The matrices are symmetric except in the
#' diagonal of Pearson, which is exactly 1 by convention. Method names
#' are propagated to \code{dimnames}.
#'
#' @keywords internal
#' @noRd
.compute_comparison_metrics <- function(cate_list) {
  m <- length(cate_list)
  nms <- names(cate_list)
  if (is.null(nms)) nms <- paste0("m", seq_len(m))
  M <- do.call(cbind, cate_list)
  rmse <- matrix(0, nrow = m, ncol = m, dimnames = list(nms, nms))
  mad  <- matrix(0, nrow = m, ncol = m, dimnames = list(nms, nms))
  pearson <- matrix(1, nrow = m, ncol = m, dimnames = list(nms, nms))
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next
      di <- M[, i] - M[, j]
      rmse[i, j] <- sqrt(mean(di^2))
      mad[i, j]  <- mean(abs(di))
      if (i < j) {
        pij <- suppressWarnings(stats::cor(M[, i], M[, j]))
        pearson[i, j] <- pij
        pearson[j, i] <- pij
      }
    }
  }
  list(rmse = rmse, pearson = pearson, mad = mad)
}

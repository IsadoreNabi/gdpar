#' Causal bridge (T-learner) between two gdpar fits
#'
#' Estimate the conditional average treatment effect (CATE) from a pair
#' of independent \code{gdpar_fit} objects fitted to disjoint arms of a
#' treatment / control design. Implements the T-learner meta-learner of
#' Kuenzel et al. (2019) on the AMM-side: each arm is fitted
#' independently and the CATE is the per-observation difference of the
#' two predictive distributions, evaluated on a common evaluation set.
#'
#' This function does not modify either fit. It assumes the two fits
#' are independent posterior samples from disjoint subsets of the
#' population (treatment arm and control arm). It does not perform any
#' causal adjustment beyond what is encoded in the two AMM
#' specifications: the assumption of no-unmeasured-confounding within
#' each arm is the responsibility of the user (see Section 4 of the
#' bridge vignette).
#'
#' @param fit_treat An object of class \code{gdpar_fit} fitted to the
#'   treatment arm.
#' @param fit_ctrl An object of class \code{gdpar_fit} fitted to the
#'   control arm. Must share the family, anchor, AMM level, and
#'   covariate structure of \code{fit_treat} (see Details).
#' @param newdata Optional data frame on which the CATE is evaluated.
#'   When \code{NULL} (default), the function attempts to recover the
#'   training data of each arm by evaluating the captured \code{data}
#'   argument of each fit's call in the caller's environment. If both
#'   recoveries succeed and the two data frames share their column
#'   structure, their \code{rbind} is used; otherwise the function
#'   aborts and requests an explicit \code{newdata}.
#' @param type Character scalar selecting the scale on which the CATE
#'   is estimated. One of \code{"response"} (default; the inverse link
#'   is applied per draw, so the CATE is on the response scale),
#'   \code{"theta_i"} (the linear predictor of the individual
#'   parameter), or \code{"linear_predictor"} (synonym of
#'   \code{"theta_i"}).
#' @param level Numeric scalar in (0, 1) with the nominal credible
#'   level for the per-observation CATE intervals. Defaults to 0.95.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return An object of class \code{gdpar_causal_bridge} with
#'   components \code{cate_draws} (matrix \code{[S, n]} when both fits
#'   are scalar, or array \code{[S, n, dim]} when both fits are
#'   multivariate or K-individual), \code{cate_mean}, \code{cate_ci},
#'   \code{newdata}, \code{id_check}, \code{fits}, \code{type},
#'   \code{level}, \code{n_draws}, \code{n_obs}, \code{call},
#'   \code{warnings} (character vector recording fallback notifications
#'   such as posterior-draw trimming; empty in the happy path), and
#'   \code{meta}. The companion S3 methods \code{print} and
#'   \code{summary} are documented in
#'   \code{\link{print.gdpar_causal_bridge}} and
#'   \code{\link{summary.gdpar_causal_bridge}}.
#'
#' @details
#' \strong{Structural compatibility.} The function aborts with
#' \code{gdpar_unsupported_feature_error} when the two fits differ in
#' any of: \code{path} (both must be the Path 1 Bayesian fit), family
#' identifier (\code{family$name}, or per-slot family identifiers when
#' \code{K > 1}), \code{K}, \code{p}, AMM level (\code{amm$level} or
#' the equivalent level inferred from each slot's spec), modulating
#' basis type (\code{amm$W$type}), anchor value, or covariate column
#' structure. The function also aborts when either fit was sampled in
#' the hierarchical regime (\code{stan_data$use_groups == 1L}); the
#' T-learner bridge for grouped fits is queued for a future sub-phase
#' and would require careful treatment of the per-group anchors which
#' is outside the scope of Sub-phase 8.5.A.
#'
#' \strong{Identifiability per arm.} The constructor records the
#' identifiability report of each fit in the \code{id_check} slot.
#' (C7) anti-aliasing of Block 6.5 is not invoked because the
#' hierarchical guard above rules out the regime in which (C7)
#' applies; this is documented for the eventual extension to grouped
#' fits.
#'
#' \strong{CATE estimator.} For each posterior draw indexed by
#' \eqn{s = 1, \dots, S} and observation \eqn{i = 1, \dots, n_{\text{new}}},
#' the bridge computes
#' \deqn{\widehat{\tau}^{(s)}_i = \hat{\mu}^{(s)}_{\text{treat}}(x_i)
#'                              - \hat{\mu}^{(s)}_{\text{ctrl}}(x_i),}
#' where \eqn{\hat{\mu}^{(s)}_{\text{arm}}(x)} is the posterior
#' prediction of the chosen \code{type} at \eqn{x}, drawn from the
#' fit's predictive distribution. The marginal posterior of the CATE
#' at each \eqn{x_i} is summarized by the empirical mean and the
#' \eqn{(\alpha/2, 1-\alpha/2)} quantiles with \eqn{\alpha = 1 - }
#' \code{level}.
#'
#' \strong{Independence of draws.} The two fits are independent (they
#' were sampled from disjoint data subsets), so the joint posterior of
#' \eqn{(\theta_{\text{treat}}, \theta_{\text{ctrl}})} factorizes and
#' any pairing of marginal draws is a valid sample from the joint.
#' The function trims to \eqn{S = \min(S_{\text{treat}}, S_{\text{ctrl}})}
#' when the two fits differ in number of draws and emits a
#' \code{gdpar_diagnostic_warning}.
#'
#' \strong{Multi-dimensional and K-individual fits.} For
#' \code{p > 1} (multivariate) and \code{K > 1} (distributional
#' regression), \code{predict.gdpar_fit} returns a 3-array of shape
#' \code{[S, n, dim]}; the CATE is computed elementwise and the
#' per-coordinate or per-slot CATEs are returned as the last
#' dimension of \code{cate_draws}. For \code{type = "response"}, the
#' canonical inverse link of each coordinate or slot is applied by
#' \code{predict.gdpar_fit} before the difference is taken; the
#' resulting CATE is therefore on the natural response scale of each
#' slot, not a uniform link-transformed scale.
#'
#' @section Methodological notes:
#' The T-learner is the most direct meta-learner to map onto the gdpar
#' pipeline: each arm is one \code{gdpar_fit} and the CATE reuses the
#' posterior machinery of \code{predict.gdpar_fit}. S-learner and
#' X-learner are queued for Block 9. The T-learner is known to suffer
#' from regularization-induced bias in unbalanced samples (see
#' Kuenzel et al. 2019, Section 3.4); the bridge vignette discusses
#' the trade-off and the alternatives.
#'
#' @section Dependencies:
#' Inherits the \pkg{posterior} dependency of \code{predict.gdpar_fit}.
#'
#' @references
#' Kuenzel, S. R., Sekhon, J. S., Bickel, P. J., and Yu, B. (2019).
#' Metalearners for estimating heterogeneous treatment effects using
#' machine learning. \emph{Proceedings of the National Academy of
#' Sciences}, 116(10), 4156-4165.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   n <- 300L
#'   df <- data.frame(x1 = rnorm(2L * n))
#'   df$arm <- rep(c("treat", "ctrl"), each = n)
#'   df$y <- with(df, ifelse(arm == "treat", 0.5, 0) + 0.8 * x1 + rnorm(2L * n, sd = 0.5))
#'   df_treat <- subset(df, arm == "treat")
#'   df_ctrl  <- subset(df, arm == "ctrl")
#'   fit_t <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df_treat,
#'                  iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   fit_c <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df_ctrl,
#'                  iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   bridge <- gdpar_causal_bridge(fit_t, fit_c,
#'                                  newdata = data.frame(x1 = seq(-2, 2, length.out = 21L)))
#'   print(bridge)
#'   summary(bridge)
#' }
#' }
#'
#' @seealso \code{\link{gdpar}}, \code{\link{predict.gdpar_fit}}
#' @export
gdpar_causal_bridge <- function(fit_treat, fit_ctrl, newdata = NULL,
                                type = c("response", "theta_i",
                                         "linear_predictor"),
                                level = 0.95, ...) {
  assert_inherits(fit_treat, "gdpar_fit", "fit_treat")
  assert_inherits(fit_ctrl, "gdpar_fit", "fit_ctrl")
  type <- match.arg(type)
  assert_numeric_scalar(level, "level", lower = 1e-3, upper = 1 - 1e-3)
  call <- match.call()
  .check_bridge_path(fit_treat, fit_ctrl)
  .check_bridge_hierarchical(fit_treat, fit_ctrl)
  .check_bridge_family(fit_treat, fit_ctrl)
  .check_bridge_dim(fit_treat, fit_ctrl)
  .check_bridge_amm(fit_treat, fit_ctrl)
  .check_bridge_anchor(fit_treat, fit_ctrl)
  newdata_resolved <- .resolve_bridge_newdata(fit_treat, fit_ctrl,
                                               newdata,
                                               eval_env = parent.frame())
  pred_t <- stats::predict(fit_treat, newdata = newdata_resolved,
                           type = type, summary = "draws")
  pred_c <- stats::predict(fit_ctrl, newdata = newdata_resolved,
                           type = type, summary = "draws")
  aligned <- .align_bridge_draws(pred_t, pred_c)
  cate_draws <- aligned$treat - aligned$ctrl
  S <- aligned$S
  alpha <- 1 - level
  ql <- alpha / 2
  qu <- 1 - alpha / 2
  summaries <- .summarize_cate(cate_draws, ql = ql, qu = qu)
  warnings_msg <- if (is.na(aligned$warning)) character(0L)
                  else aligned$warning
  obj <- list(
    cate_draws = cate_draws,
    cate_mean  = summaries$mean,
    cate_ci    = summaries$ci,
    newdata    = newdata_resolved,
    id_check   = list(
      treat = fit_treat$identifiability_report,
      ctrl  = fit_ctrl$identifiability_report
    ),
    fits       = list(treat = fit_treat, ctrl = fit_ctrl),
    type       = type,
    level      = level,
    n_draws    = S,
    n_obs      = nrow(newdata_resolved),
    call       = call,
    warnings   = warnings_msg,
    meta       = list(
      dim_kind = summaries$dim_kind,
      dim_size = summaries$dim_size,
      dim_names = summaries$dim_names,
      newdata_source = attr(newdata_resolved, "bridge_source")
    )
  )
  class(obj) <- c("gdpar_causal_bridge", "list")
  obj
}

#' Internal: assert both fits use the Path 1 Bayesian fit
#' @keywords internal
#' @noRd
.check_bridge_path <- function(fit_treat, fit_ctrl) {
  path_t <- fit_treat$path %||% "bayes"
  path_c <- fit_ctrl$path %||% "bayes"
  if (!identical(path_t, "bayes") || !identical(path_c, "bayes")) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to be Path 1",
          "(path = 'bayes'); received fit_treat$path = '%s' and",
          "fit_ctrl$path = '%s'. T-learner support for Paths 2/3 is",
          "queued for a future sub-phase."
        ),
        path_t, path_c
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(path_treat = path_t, path_ctrl = path_c)
    )
  }
  invisible(NULL)
}

#' Internal: assert neither fit is hierarchical (use_groups == 0 on both)
#'
#' Mirrors the pattern of gdpar_bvm_check (Session 1 of Phase E,
#' decision 3.4): the T-learner bridge for fits with a per-group
#' anchor introduces a random-effect dimension whose handling under
#' the T-learner difference operator is not standard and is queued for
#' a later sub-phase. The abort uses the same condition class used by
#' bvm_check so user code that already handles unsupported-feature
#' errors covers both helpers.
#'
#' @keywords internal
#' @noRd
.check_bridge_hierarchical <- function(fit_treat, fit_ctrl) {
  is_grouped <- function(fit) {
    !is.null(fit$stan_data$use_groups) &&
      as.integer(fit$stan_data$use_groups) == 1L
  }
  if (is_grouped(fit_treat) || is_grouped(fit_ctrl)) {
    gdpar_abort(
      paste(
        "gdpar_causal_bridge does not apply to hierarchical fits",
        "(Block 6.5 grouped anchor regime). The T-learner difference of",
        "per-group anchors is not defined in the canonical formulation",
        "of the bridge; per-group bridges are queued for a future",
        "sub-phase. Refit each arm without the 'group' argument to",
        "build a bridge on the non-hierarchical regime."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(
        treat_grouped = is_grouped(fit_treat),
        ctrl_grouped  = is_grouped(fit_ctrl)
      )
    )
  }
  invisible(NULL)
}

#' Internal: assert the two fits have compatible family identifiers
#' @keywords internal
#' @noRd
.check_bridge_family <- function(fit_treat, fit_ctrl) {
  fam_t <- fit_treat$family
  fam_c <- fit_ctrl$family
  name_t <- fam_t$name
  name_c <- fam_c$name
  link_t <- fam_t$link
  link_c <- fam_c$link
  if (!identical(name_t, name_c) || !identical(link_t, link_c)) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to share the family",
          "and link; received treat = %s (link = %s) and",
          "ctrl = %s (link = %s)."
        ),
        name_t, link_t, name_c, link_c
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(
        family_treat = name_t, family_ctrl = name_c,
        link_treat = link_t, link_ctrl = link_c
      )
    )
  }
  ps_t <- fam_t$param_specs
  ps_c <- fam_c$param_specs
  if (!is.null(ps_t) || !is.null(ps_c)) {
    n_t <- length(ps_t %||% list())
    n_c <- length(ps_c %||% list())
    if (n_t != n_c) {
      gdpar_abort(
        sprintf(
          paste(
            "gdpar_causal_bridge requires both fits to share the family",
            "structure; received %d slot(s) in treat and %d in ctrl."
          ),
          n_t, n_c
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(n_slots_treat = n_t, n_slots_ctrl = n_c)
      )
    }
    names_t <- vapply(ps_t, function(s) as.character(s$family_id %||% s$name %||% NA_character_),
                      character(1L))
    names_c <- vapply(ps_c, function(s) as.character(s$family_id %||% s$name %||% NA_character_),
                      character(1L))
    if (!identical(names_t, names_c)) {
      gdpar_abort(
        sprintf(
          paste(
            "gdpar_causal_bridge requires both fits to share the per-slot",
            "family identifiers; received treat = (%s) and ctrl = (%s)."
          ),
          paste(names_t, collapse = ", "),
          paste(names_c, collapse = ", ")
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(slot_families_treat = names_t,
                    slot_families_ctrl  = names_c)
      )
    }
  }
  invisible(NULL)
}

#' Internal: assert the two fits have the same dim (K and p)
#' @keywords internal
#' @noRd
.check_bridge_dim <- function(fit_treat, fit_ctrl) {
  K_t <- fit_treat[["K"]] %||% 1L
  K_c <- fit_ctrl[["K"]] %||% 1L
  p_t <- fit_treat[["p"]] %||% 1L
  p_c <- fit_ctrl[["p"]] %||% 1L
  if (K_t != K_c || p_t != p_c) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to have matching",
          "(K, p); received treat (K = %d, p = %d) and ctrl",
          "(K = %d, p = %d)."
        ),
        K_t, p_t, K_c, p_c
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(K_treat = K_t, K_ctrl = K_c,
                  p_treat = p_t, p_ctrl = p_c)
    )
  }
  invisible(NULL)
}

#' Internal: assert AMM compatibility (level and modulating basis type)
#'
#' Compares the AMM level and the modulating basis type. The level
#' compares the structural composition of the spec (a / b / W); the
#' basis type ensures the predict path on newdata reuses the same
#' algorithm (polynomial vs B-spline).
#'
#' @keywords internal
#' @noRd
.check_bridge_amm <- function(fit_treat, fit_ctrl) {
  amm_t <- fit_treat$amm %||% fit_treat$amm_list_canonical
  amm_c <- fit_ctrl$amm %||% fit_ctrl$amm_list_canonical
  if (is.null(amm_t) || is.null(amm_c)) {
    gdpar_abort(
      "Internal error: AMM spec missing from one of the fits.",
      class = "gdpar_internal_error"
    )
  }
  level_t <- .bridge_amm_level(amm_t)
  level_c <- .bridge_amm_level(amm_c)
  if (!identical(level_t, level_c)) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to share the AMM",
          "level (or per-slot levels when K > 1); received treat = %s",
          "and ctrl = %s."
        ),
        paste(level_t, collapse = ", "),
        paste(level_c, collapse = ", ")
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(level_treat = level_t, level_ctrl = level_c)
    )
  }
  W_type_t <- .bridge_amm_W_type(amm_t)
  W_type_c <- .bridge_amm_W_type(amm_c)
  if (!identical(W_type_t, W_type_c)) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to share the",
          "modulating basis type; received treat = %s and ctrl = %s."
        ),
        W_type_t %||% "<none>", W_type_c %||% "<none>"
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(W_type_treat = W_type_t, W_type_ctrl = W_type_c)
    )
  }
  cov_t <- .bridge_amm_covariates(amm_t)
  cov_c <- .bridge_amm_covariates(amm_c)
  if (!identical(cov_t, cov_c)) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to share the",
          "covariate column structure of the AMM design (a, b, W).",
          "Mismatched components: %s."
        ),
        paste(setdiff(union(names(cov_t), names(cov_c)),
                      intersect(names(cov_t), names(cov_c))),
              collapse = ", ")
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(covariates_treat = cov_t, covariates_ctrl = cov_c)
    )
  }
  invisible(NULL)
}

#' Internal: assert numerically close anchors
#'
#' The anchor enters the modulating term as theta_ref^k - anchor^k,
#' so a mismatch between arms changes the meaning of the predicted
#' theta_i and therefore of the CATE. The check tolerates small
#' numerical noise (relative or absolute tolerance of 1e-8).
#'
#' @keywords internal
#' @noRd
.check_bridge_anchor <- function(fit_treat, fit_ctrl) {
  a_t <- fit_treat$anchor
  a_c <- fit_ctrl$anchor
  if (length(a_t) != length(a_c)) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to share the anchor",
          "length; received treat (length = %d) and ctrl (length = %d)."
        ),
        length(a_t), length(a_c)
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(anchor_treat_length = length(a_t),
                  anchor_ctrl_length  = length(a_c))
    )
  }
  diff <- abs(a_t - a_c)
  scale <- pmax(abs(a_t), abs(a_c), 1)
  if (any(diff > 1e-8 * scale)) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge requires both fits to share the anchor",
          "value; received treat = (%s) and ctrl = (%s). Refit one arm",
          "anchored to the other arm's value."
        ),
        paste(format(a_t, digits = 6), collapse = ", "),
        paste(format(a_c, digits = 6), collapse = ", ")
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(anchor_treat = a_t, anchor_ctrl = a_c)
    )
  }
  invisible(NULL)
}

#' Internal: infer the AMM level of an amm_spec or list of amm_spec
#' @keywords internal
#' @noRd
.bridge_amm_level <- function(amm) {
  if (inherits(amm, "amm_spec")) {
    return(as.integer(amm$level))
  }
  if (is.list(amm)) {
    return(vapply(amm, function(a) {
      if (inherits(a, "amm_spec")) as.integer(a$level) else NA_integer_
    }, integer(1L)))
  }
  NA_integer_
}

#' Internal: extract the W basis type ("polynomial" / "bspline") or NULL
#' @keywords internal
#' @noRd
.bridge_amm_W_type <- function(amm) {
  if (inherits(amm, "amm_spec")) {
    if (is.null(amm$W)) return(NULL)
    return(as.character(amm$W$type))
  }
  if (is.list(amm)) {
    types <- vapply(amm, function(a) {
      if (inherits(a, "amm_spec") && !is.null(a$W)) {
        as.character(a$W$type)
      } else {
        NA_character_
      }
    }, character(1L))
    if (all(is.na(types))) return(NULL)
    return(types)
  }
  NULL
}

#' Internal: extract the covariate names per AMM component
#' @keywords internal
#' @noRd
.bridge_amm_covariates <- function(amm) {
  one <- function(a) {
    list(
      a_vars = if (!is.null(a$a)) all.vars(a$a) else character(0L),
      b_vars = if (!is.null(a$b)) all.vars(a$b) else character(0L),
      x_vars = if (!is.null(a$x_vars)) a$x_vars else character(0L)
    )
  }
  if (inherits(amm, "amm_spec")) {
    return(one(amm))
  }
  if (is.list(amm)) {
    return(lapply(amm, one))
  }
  list()
}

#' Internal: resolve newdata for the bridge
#'
#' Returns the supplied newdata when non-NULL, marking
#' \code{attr(., "bridge_source") = "user"}. When NULL, attempts to
#' recover both arms' training data by evaluating the captured
#' \code{data} argument of each fit's call in \code{eval_env}. If both
#' recoveries succeed and the two data frames share their column
#' names, the function returns their \code{rbind} marked with
#' \code{attr(., "bridge_source") = "training_rbind"}; otherwise it
#' aborts with a class-tagged error and instructs the user to pass
#' \code{newdata} explicitly.
#'
#' @keywords internal
#' @noRd
.resolve_bridge_newdata <- function(fit_treat, fit_ctrl, newdata,
                                     eval_env) {
  if (!is.null(newdata)) {
    assert_data_frame(newdata, "newdata")
    attr(newdata, "bridge_source") <- "user"
    return(newdata)
  }
  recover <- function(fit, arm_label) {
    call <- fit$call
    if (is.null(call) || is.null(call$data)) {
      return(NULL)
    }
    tryCatch(
      eval(call$data, envir = eval_env),
      error = function(e) NULL
    )
  }
  data_t <- recover(fit_treat, "treat")
  data_c <- recover(fit_ctrl, "ctrl")
  if (is.null(data_t) || is.null(data_c) ||
      !is.data.frame(data_t) || !is.data.frame(data_c)) {
    gdpar_abort(
      paste(
        "gdpar_causal_bridge cannot recover the training data of one",
        "or both fits from the captured call. Pass an explicit",
        "'newdata' argument with the evaluation grid for the CATE."
      ),
      class = "gdpar_input_error",
      data = list(treat_recovered = !is.null(data_t),
                  ctrl_recovered  = !is.null(data_c))
    )
  }
  if (!identical(sort(colnames(data_t)), sort(colnames(data_c)))) {
    gdpar_abort(
      sprintf(
        paste(
          "gdpar_causal_bridge recovered training data for both arms",
          "but their column structures do not match (treat = (%s),",
          "ctrl = (%s)). Pass an explicit 'newdata' argument."
        ),
        paste(colnames(data_t), collapse = ", "),
        paste(colnames(data_c), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(colnames_treat = colnames(data_t),
                  colnames_ctrl  = colnames(data_c))
    )
  }
  out <- rbind(data_t[, colnames(data_t), drop = FALSE],
               data_c[, colnames(data_t), drop = FALSE])
  attr(out, "bridge_source") <- "training_rbind"
  out
}

#' Internal: align the per-draw arrays of the two arms
#'
#' Trims the longer of the two draws sets to the length of the shorter,
#' emitting a \code{gdpar_diagnostic_warning} when trimming occurs.
#' Returns a list with components \code{treat}, \code{ctrl}, \code{S},
#' and \code{warning} (character scalar with the trimming notification
#' when trimming occurred, \code{NA_character_} otherwise). The two
#' inputs must be either both matrices or both arrays of the same
#' shape (modulo the first axis). The \code{warning} field is persisted
#' by the constructor into the \code{$warnings} slot of the resulting
#' \code{gdpar_causal_bridge} object so the print method can surface
#' the fallback notification (D48 canonical norm; D50 of Sesion 18
#' Etapa 2 of Sesion 8.4).
#'
#' @keywords internal
#' @noRd
.align_bridge_draws <- function(pred_t, pred_c) {
  S_t <- if (is.matrix(pred_t)) nrow(pred_t) else dim(pred_t)[1L]
  S_c <- if (is.matrix(pred_c)) nrow(pred_c) else dim(pred_c)[1L]
  S <- min(S_t, S_c)
  warning_msg <- NA_character_
  if (S_t != S_c) {
    warning_msg <- sprintf(
      paste(
        "gdpar_causal_bridge trimmed posterior draws to S = %d",
        "(treat had %d, ctrl had %d)."
      ),
      S, S_t, S_c
    )
    gdpar_warn(
      warning_msg,
      class = "gdpar_diagnostic_warning",
      data = list(S_treat = S_t, S_ctrl = S_c, S = S)
    )
  }
  trim_first_axis <- function(arr, S) {
    if (is.matrix(arr)) {
      return(arr[seq_len(S), , drop = FALSE])
    }
    nd <- length(dim(arr))
    idx <- c(list(seq_len(S)), rep(list(quote(expr =)), nd - 1L))
    do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
  }
  list(
    treat   = trim_first_axis(pred_t, S),
    ctrl    = trim_first_axis(pred_c, S),
    S       = S,
    warning = warning_msg
  )
}

#' Internal: summarize the CATE draws into mean and credible intervals
#'
#' Returns a list with the per-observation posterior mean and the
#' \code{(ql, qu)} quantiles. The shape of the outputs follows the
#' shape of the input draws: matrix \code{[S, n]} produces a vector of
#' length n and a matrix \code{[n, 2]}; array \code{[S, n, dim]}
#' produces a matrix \code{[n, dim]} and an array \code{[n, dim, 2]}.
#' The metadata fields \code{dim_kind}, \code{dim_size},
#' \code{dim_names} record the structure for downstream methods.
#'
#' @keywords internal
#' @noRd
.summarize_cate <- function(cate_draws, ql, qu) {
  if (is.matrix(cate_draws)) {
    n <- ncol(cate_draws)
    cate_mean <- colMeans(cate_draws)
    q <- apply(cate_draws, 2L, stats::quantile,
                probs = c(ql, qu), names = FALSE)
    cate_ci <- matrix(NA_real_, nrow = n, ncol = 2L)
    cate_ci[, 1L] <- q[1L, ]
    cate_ci[, 2L] <- q[2L, ]
    colnames(cate_ci) <- c("lower", "upper")
    return(list(
      mean = cate_mean,
      ci   = cate_ci,
      dim_kind = "scalar",
      dim_size = 1L,
      dim_names = NULL
    ))
  }
  if (length(dim(cate_draws)) == 3L) {
    n <- dim(cate_draws)[2L]
    K <- dim(cate_draws)[3L]
    slot_names <- dimnames(cate_draws)[[3L]]
    cate_mean <- apply(cate_draws, c(2L, 3L), mean)
    cate_ci_lower <- apply(cate_draws, c(2L, 3L), stats::quantile,
                            probs = ql, names = FALSE)
    cate_ci_upper <- apply(cate_draws, c(2L, 3L), stats::quantile,
                            probs = qu, names = FALSE)
    cate_ci <- array(NA_real_, dim = c(n, K, 2L),
                      dimnames = list(NULL, slot_names,
                                       c("lower", "upper")))
    cate_ci[, , 1L] <- cate_ci_lower
    cate_ci[, , 2L] <- cate_ci_upper
    dimnames(cate_mean) <- list(NULL, slot_names)
    dim_kind <- if (is.null(slot_names) ||
                    all(grepl("^dim_", slot_names))) {
      "multi"
    } else {
      "K_individual"
    }
    return(list(
      mean = cate_mean,
      ci   = cate_ci,
      dim_kind = dim_kind,
      dim_size = K,
      dim_names = slot_names
    ))
  }
  gdpar_abort(
    "Internal error: unsupported shape for cate_draws.",
    class = "gdpar_internal_error",
    data = list(dim = dim(cate_draws))
  )
}


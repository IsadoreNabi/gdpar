#' Internal: deduplicate consecutive message blocks captured from stderr
#'
#' Splits a character vector of captured stderr lines into blocks
#' separated by empty lines, and returns a flattened character vector
#' that preserves the first occurrence of each unique block while
#' dropping subsequent identical blocks. Used to silence the duplicate
#' divergence warnings emitted by cmdstanr when \code{verbose = FALSE}.
#'
#' @param lines Character vector of captured stderr lines.
#' @return Character vector ready to be re-emitted via \code{message()}.
#' @keywords internal
#' @noRd
dedup_message_blocks <- function(lines) {
  if (length(lines) == 0L) return(character(0))
  blocks <- list()
  current <- character(0)
  for (line in lines) {
    if (nzchar(line)) {
      current <- c(current, line)
    } else if (length(current) > 0L) {
      blocks[[length(blocks) + 1L]] <- current
      current <- character(0)
    }
  }
  if (length(current) > 0L) {
    blocks[[length(blocks) + 1L]] <- current
  }
  seen_keys <- character(0)
  out <- character(0)
  for (b in blocks) {
    key <- paste(b, collapse = "\n")
    if (!key %in% seen_keys) {
      seen_keys <- c(seen_keys, key)
      out <- c(out, b, "")
    }
  }
  if (length(out) > 0L && !nzchar(out[length(out)])) {
    out <- out[-length(out)]
  }
  out
}

#' Fit an AMM canonical model
#'
#' Main entry point of the package. Fits the AMM canonical
#' decomposition of the individual parameter via one of three
#' estimation paths. In this version of the package, only Path 1
#' (hierarchical Bayesian inference via Stan) is fully implemented;
#' Path 2 and Path 3 wrappers are placeholders that signal the
#' development status.
#'
#' @param formula Either a two-sided formula or an object of class
#'   \code{\link{gdpar_formula_set}}. In the legacy two-sided form,
#'   the left-hand side is the outcome variable; the right-hand side
#'   either lists the covariates that enter the modulating component
#'   as the linear factor x in \eqn{W(\theta) x} (the AMM components
#'   are declared via the \code{amm} argument), or it contains the
#'   AMM wrapper calls \code{a(...)}, \code{b(...)} and/or \code{W()}
#'   themselves (canonized in sub-phase 8.3.3 as the K = 1 entry to
#'   the high-level parser path). A \code{gdpar_formula_set}
#'   (constructed via \code{\link{gdpar_formula_set}} or via the
#'   brms-style sugar \code{\link{gdpar_bf}}) is the multi-parameter
#'   form: one slot per K-individual parameter, with AMM wrappers in
#'   each slot's right-hand side.
#' @param W An object of class \code{W_basis} produced by
#'   \code{\link{W_basis}}, or \code{NULL} (default). Used by the
#'   K-individual paths (formula set, named list of \code{amm_spec},
#'   or classic formula with AMM wrappers in the RHS) to supply the
#'   modulating basis when one or more slots declare \code{W()}.
#'   Ignored in the legacy single-\code{amm_spec} path, where the
#'   modulating basis travels through \code{amm$W}.
#' @param family An object of class \code{gdpar_family} produced by
#'   \code{\link{gdpar_family}} or \code{\link{gdpar_family_custom}}.
#'   Defaults to \code{gdpar_family("gaussian")}.
#' @param amm Either an object of class \code{amm_spec} produced by
#'   \code{\link{amm_spec}} (legacy single-parameter path, K = 1) or
#'   a named list of \code{amm_spec} objects (K-individual low-level
#'   path closed in sub-phase 8.3.3, names matching the eligible
#'   parameters of \code{family}). Defaults to a Level 0 (degenerate)
#'   single \code{amm_spec} when not supplied. When \code{formula} is
#'   a \code{\link{gdpar_formula_set}} or contains AMM wrapper calls
#'   (\code{a()}/\code{b()}/\code{W()}) in its right-hand side, the
#'   \code{amm} argument must remain at its default; the canonical
#'   specifications come from the formula in those paths.
#' @param data A data frame containing the variables referenced by
#'   \code{formula} and \code{amm}.
#' @param prior An object of class \code{gdpar_prior} produced by
#'   \code{\link{gdpar_prior}}. When NULL (default), the package
#'   defaults are used.
#' @param path Character scalar identifying the estimation path. One
#'   of \code{"bayes"} (Path 1, default), \code{"vcm"} (Path 2),
#'   \code{"hyper"} (Path 3).
#' @param anchor Either a numeric scalar with the anchor value
#'   theta_0 in the linear-predictor scale, or one of
#'   \code{"prior_mean"} (default) and \code{"empirical_y"} (the
#'   linkfun applied to the outcome mean).
#' @param skip_id_check Logical scalar. If \code{TRUE}, skip the
#'   automatic Gram-matrix identifiability diagnostic. Defaults to
#'   \code{FALSE}; turning it off is intended for advanced users who
#'   have verified identifiability separately.
#' @param chains Integer scalar with the number of Hamiltonian Monte
#'   Carlo chains. Defaults to 4.
#' @param iter_warmup Integer scalar with the number of warmup
#'   iterations per chain. Defaults to 1000.
#' @param iter_sampling Integer scalar with the number of sampling
#'   iterations per chain. Defaults to 1000.
#' @param adapt_delta Numeric scalar in (0, 1) controlling the target
#'   acceptance probability of the No-U-Turn Sampler. Defaults to 0.95
#'   (more conservative than the Stan default 0.8) to reduce divergent
#'   transitions in hierarchical models.
#' @param max_treedepth Integer scalar bounding the tree depth of the
#'   No-U-Turn Sampler. Defaults to 12.
#' @param refresh Integer scalar with the progress-reporting frequency
#'   passed to \code{cmdstanr}. Defaults to 100, which prints a status
#'   line every 100 iterations.
#' @param verbose Logical scalar controlling the verbosity of
#'   informational messages from the package itself (independent of
#'   \code{refresh}). Defaults to \code{TRUE}.
#' @param seed Optional integer scalar passed to \code{cmdstanr} for
#'   reproducibility. Defaults to \code{NULL}, which lets
#'   \code{cmdstanr} pick a seed. The same seed is forwarded to the
#'   pre-flight diagnostic when \code{parametrization = "auto"}.
#' @param group Optional one-sided formula identifying the grouping
#'   variable in \code{data} that anchors per-group hierarchical
#'   estimation of \code{theta_ref}. When \code{NULL} (default) the
#'   model uses a single global \code{theta_ref}, matching the Block 6
#'   behavior bit-exactly. When supplied, \code{theta_ref} is promoted
#'   to a vector indexed by group with a Normal hyperprior
#'   \eqn{\theta_{ref}[g] \sim \mathrm{Normal}(\mu_{\theta_{ref}},
#'   \sigma_{\theta_{ref}})}, analogous to a random intercept. Only
#'   one-sided formulas with a single variable name are accepted (for
#'   example \code{~ species}); the variable must exist in \code{data}.
#'   The grouping is rejected at pre-flight if (C7) of Block 6.5 is
#'   violated: that is, if \code{a} or \code{b} contain columns that are
#'   constant within every group level or otherwise rank-deficient with
#'   the group indicator (perfect aliasing with the per-group anchor).
#' @param parametrization Character scalar selecting the sampling
#'   parametrization for the additive component \code{a} and the
#'   modulating component \code{W}. One of \code{"auto"} (default,
#'   runs a short pre-flight NCP fit and decides per-component via a
#'   sequence of three filters: a divergence-attribution t-test, an
#'   E-BFMI threshold, and a chain-aware block-bootstrap z-test on
#'   the posterior-to-prior contraction of the effective coefficient,
#'   Path B'), \code{"ncp"} (forces non-centered parametrization for
#'   both components, skipping the pre-flight), \code{"cp"} (forces
#'   centered parametrization for both components, skipping the
#'   pre-flight). The pre-flight adds approximately 30% wall-time per
#'   fit when active. See
#'   \code{vignette("vop01_parametrization_toggle", package = "gdpar")}
#'   for the operational guide.
#' @param parametrization_a Optional character scalar
#'   (\code{"ncp"} or \code{"cp"}) overriding the parametrization for
#'   component \code{a}. When \code{NULL} (default), inherits from
#'   \code{parametrization}.
#' @param parametrization_aggregation Optional character scalar
#'   (\code{"any_ncp"}, \code{"majority"}, or \code{"per_k"})
#'   controlling how the per-coordinate CP/NCP decisions of the
#'   multivariate (\code{amm$p > 1L}) pre-flight are aggregated to a
#'   per-component decision. \code{"any_ncp"} (default) is
#'   conservative: a component is CP only when every coordinate's
#'   decision is CP. \code{"majority"} selects CP when the strict
#'   majority of coordinates votes CP; ties break toward NCP.
#'   \code{"per_k"} keeps the per-coordinate decisions intact and
#'   honors them at the sampling level via per-k
#'   \code{segment()}-based priors (Phase H.2 of Block 5.2). The
#'   uniform branches (all-CP or all-NCP) compile to byte-identical
#'   Stan code with the H.1 template, preserving bit-exact paridad.
#'   Ignored for the univariate path (\code{amm$p == 1L}). The full
#'   per-coord report is stored in \code{fit$parametrization$report}
#'   as an object of class \code{gdpar_preflight_report}.
#' @param parametrization_W Optional character scalar
#'   (\code{"ncp"} or \code{"cp"}) overriding the parametrization for
#'   component \code{W}. When \code{NULL} (default), inherits from
#'   \code{parametrization}. If both \code{parametrization_a} and
#'   \code{parametrization_W} are explicit, the pre-flight is
#'   skipped regardless of \code{parametrization}.
#' @param id_check_rigor Character scalar, one of \code{"full"}
#'   (default) or \code{"fast"}, controlling the strictness of the
#'   basis-restricted identifiability check. \code{"full"} aborts on
#'   any C1-C4 (univariate path) or C4-bis cross-coordinate
#'   identifiability violation (multivariate path). \code{"fast"}
#'   downgrades the C4-bis check to a single consolidated warning at
#'   the end of the per-coordinate loop, allowing the fit to proceed.
#'   The design pattern for legitimate overlap between the additive
#'   and modulating channels is documented in
#'   \code{vignette("v01_amm_identifiability", package = "gdpar")}
#'   Section 6.6.1.7. \code{"fast"} is only meaningful for the
#'   multivariate path (\code{amm$p > 1L}); for the univariate path
#'   the value is forwarded but the C4-bis check is not run.
#' @param ... Additional arguments forwarded to the underlying
#'   sampler.
#'
#' @return An object of class \code{gdpar_fit} with components
#'   \code{fit} (the underlying \code{cmdstanr} fit object),
#'   \code{amm}, \code{family}, \code{prior}, \code{design},
#'   \code{anchor}, \code{identifiability_report}, \code{diagnostics},
#'   \code{parametrization} (a list with the resolved CP/NCP flags
#'   per component plus the pre-flight diagnostic statistics when
#'   applicable), \code{call} and \code{path}.
#'
#' @details
#' The function orchestrates the Path 1 fit in five steps: (1) input
#' validation and standardization of the covariates entering the
#' design matrices; (2) the basis-restricted identifiability
#' diagnostic of \code{\link{gdpar_check_identifiability}}; (3)
#' generation of the Stan model source by substituting the prior
#' placeholders into the static template
#' \code{inst/stan/amm_main.stan}; (4) compilation and sampling via
#' \code{cmdstanr}; (5) collection of convergence diagnostics from
#' \pkg{posterior} and assembly of the returned object.
#'
#' Covariates entering the additive and multiplicative bases are
#' centered before fitting, enforcing assumption (C1) of Block 1 at
#' the empirical level. The covariates entering the modulating
#' component are additionally scaled to unit standard deviation; this
#' standardization is recorded so that predictions on new data can
#' apply the same transformation consistently.
#'
#' The anchor value theta_0 enters the parametrization
#' W(theta) - W(theta_anchor) and is therefore a parametrization
#' device, not an inferential statement about the posterior of
#' theta_ref. Choosing a different anchor changes the parametrization
#' but not the data-generating model.
#'
#' @section Methodological notes:
#' Path 1 supports finite-dimensional parametric AMM specifications
#' in this version. Non-parametric extensions (Gaussian-process
#' priors on a, b, W; adaptive splines with growing basis dimension)
#' are deferred to a future version; users who need such flexibility
#' are referred to Path 2 (varying-coefficient models via mgcv).
#'
#' The polynomial restriction on \code{W_basis} in Path 1 is also a
#' v0 limitation; B-spline and user-defined W bases require either a
#' Stan-side basis evaluator that the current template does not
#' implement, or a precomputed basis values approach that complicates
#' the data block. Both extensions are planned.
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} (Suggests) for compilation and sampling, and
#' \pkg{posterior} (Imports) for convergence diagnostics on the
#' resulting posterior draws. \pkg{cmdstanr} is loaded conditionally
#' via \code{requireNamespace}; the function aborts with an
#' informative error if the package is not installed.
#'
#' @references
#' See \code{vignette("v01_amm_identifiability", package = "gdpar")}
#' for the canonical form and identifiability conditions;
#' \code{vignette("v04_asymptotics_path1_bayesian", package = "gdpar")}
#' for the asymptotic theory of Path 1; and
#' \code{vignette("vop01_parametrization_toggle", package = "gdpar")}
#' for the operational guide to the CP/NCP toggle, including the
#' three-filter pre-flight diagnostic, introspection of
#' \code{fit$parametrization$meta}, and known limitations under
#' confounding.
#'
#' Stan Development Team (2024). Stan User's Guide, version 2.35.
#'
#' Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
#' Betancourt, M., Brubaker, M., Guo, J., Li, P., and Riddell, A.
#' (2017). Stan: A probabilistic programming language. Journal of
#' Statistical Software, 76(1).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   set.seed(NULL)
#'   n <- 200
#'   df <- data.frame(
#'     x1 = rnorm(n),
#'     x2 = rnorm(n)
#'   )
#'   df$y <- with(df, 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n, sd = 0.5))
#'   spec <- amm_spec(a = ~ x1 + x2)
#'   fit <- gdpar(
#'     formula       = y ~ x1 + x2,
#'     family        = gdpar_family("gaussian"),
#'     amm           = spec,
#'     data          = df,
#'     iter_warmup   = 200,
#'     iter_sampling = 200,
#'     chains        = 2
#'   )
#'   print(fit)
#' }
#' }
#'
#' @seealso \code{\link{amm_spec}}, \code{\link{W_basis}},
#'   \code{\link{gdpar_family}}, \code{\link{gdpar_prior}},
#'   \code{\link{gdpar_check_identifiability}},
#'   \code{\link{gdpar_bvm_check}},
#'   \code{\link{gdpar_contraction_diagnostic}}
#' @export
gdpar <- function(formula, family = gdpar_family("gaussian"),
                  amm = amm_spec(),
                  W = NULL,
                  data,
                  prior = NULL,
                  path = c("bayes", "vcm", "hyper"),
                  anchor = "prior_mean",
                  skip_id_check = FALSE,
                  chains = 4L,
                  iter_warmup = 1000L,
                  iter_sampling = 1000L,
                  adapt_delta = 0.95,
                  max_treedepth = 12L,
                  refresh = 100L,
                  verbose = TRUE,
                  seed = NULL,
                  group = NULL,
                  parametrization = c("auto", "ncp", "cp"),
                  parametrization_a = NULL,
                  parametrization_W = NULL,
                  parametrization_aggregation = NULL,
                  id_check_rigor = c("full", "fast"), ...) {
  path <- match.arg(path)
  parametrization <- match.arg(parametrization)
  id_check_rigor <- match.arg(id_check_rigor)
  if (!is.null(parametrization_a)) {
    parametrization_a <- match.arg(parametrization_a,
                                    c("ncp", "cp"))
  }
  if (!is.null(parametrization_W)) {
    parametrization_W <- match.arg(parametrization_W,
                                    c("ncp", "cp"))
  }
  if (!is.null(parametrization_aggregation)) {
    parametrization_aggregation <- match.arg(
      parametrization_aggregation,
      c("any_ncp", "majority", "per_k")
    )
  }
  call <- match.call()
  if (path == "vcm") {
    gdpar_abort(
      paste(
        "Path 2 (varying-coefficient via mgcv) is planned for a future",
        "version. Use path = 'bayes' for the hierarchical Bayesian fit",
        "or fit a varying-coefficient model directly with",
        "mgcv::gam(...)."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  if (path == "hyper") {
    gdpar_abort(
      paste(
        "Path 3 (hypernetwork via torch) is planned for a future",
        "release; torch is not currently a declared dependency of the",
        "package. Use path = 'bayes' for the hierarchical Bayesian fit."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  .family_is_named_list <- is.list(family) &&
                           !inherits(family, "gdpar_family") &&
                           !inherits(family, "gdpar_family_multi") &&
                           !is.null(names(family)) &&
                           all(nzchar(names(family))) &&
                           anyDuplicated(names(family)) == 0L &&
                           all(vapply(family, inherits,
                                      logical(1L), "gdpar_family"))
  if (!.family_is_named_list) {
    assert_inherits(family, c("gdpar_family", "gdpar_family_multi"),
                    "family")
  }
  .family_input_raw <- family
  .formula_set_input <- inherits(formula, "gdpar_formula_set")
  .amm_list_input <- is.list(amm) && !inherits(amm, "amm_spec") &&
                      !is.null(names(amm))
  .classic_with_amm_calls <- !.formula_set_input &&
                              inherits(formula, "formula") &&
                              length(formula) == 3L &&
                              .gdpar_rhs_has_amm_calls(formula)
  if (.formula_set_input || .amm_list_input || .classic_with_amm_calls) {
    if (.formula_set_input) {
      if (!.gdpar_is_default_amm_spec(amm)) {
        gdpar_abort(
          paste(
            "When 'formula' is a gdpar_formula_set, the 'amm' argument",
            "must remain at its default amm_spec(): the formula set is the",
            "canonical declaration of every AMM component per slot."
          ),
          class = "gdpar_input_error"
        )
      }
      amm_list_canonical <- .gdpar_formula_set_to_amm_spec_list(formula, W)
      outcome_name <- formula$outcome
      formula_env <- formula$env
    } else if (.amm_list_input) {
      amm_list_canonical <- amm
      for (slot in names(amm_list_canonical)) {
        if (!nzchar(slot)) {
          gdpar_abort(
            paste(
              "Every entry of the named list passed to 'amm' must carry a",
              "non-empty slot name (the canonical parameter name)."
            ),
            class = "gdpar_input_error"
          )
        }
        if (!inherits(amm_list_canonical[[slot]], "amm_spec")) {
          gdpar_abort(
            sprintf(
              paste(
                "Entry '%s' of the named list passed to 'amm' must be an",
                "object of class amm_spec; received an object of class %s."
              ),
              slot,
              paste(sQuote(class(amm_list_canonical[[slot]])),
                    collapse = ", ")
            ),
            class = "gdpar_input_error",
            data = list(slot = slot,
                        received = class(amm_list_canonical[[slot]]))
          )
        }
      }
      if (anyDuplicated(names(amm_list_canonical)) > 0L) {
        gdpar_abort(
          "Named list passed to 'amm' has duplicated slot names.",
          class = "gdpar_input_error"
        )
      }
      if (!inherits(formula, "formula") || length(formula) != 3L) {
        gdpar_abort(
          paste(
            "When 'amm' is a named list of amm_spec, 'formula' must be a",
            "two-sided formula 'y ~ ...'."
          ),
          class = "gdpar_input_error"
        )
      }
      outcome_name <- as.character(formula[[2L]])
      formula_env <- environment(formula)
    } else {
      if (!.gdpar_is_default_amm_spec(amm)) {
        gdpar_abort(
          paste(
            "The right-hand side of 'formula' contains AMM wrapper calls",
            "(`a()`/`b()`/`W()`). In this path the 'amm' argument must",
            "remain at its default amm_spec(); the wrappers in the RHS",
            "are the canonical declaration of the AMM components."
          ),
          class = "gdpar_input_error"
        )
      }
      first_eligible <- family$param_specs[[1L]]$name
      args_for_fs <- list(formula)
      names(args_for_fs) <- first_eligible
      fs <- do.call(gdpar_formula_set, args_for_fs)
      amm_list_canonical <- .gdpar_formula_set_to_amm_spec_list(fs, W)
      outcome_name <- fs$outcome
      formula_env <- fs$env
    }
    K <- length(amm_list_canonical)
    if (K > 1L) {
      if (.family_is_named_list) {
        .het_resolved <- .gdpar_resolve_heterogeneous_family_K(
          .family_input_raw, names(amm_list_canonical)
        )
        family_promoted     <- .het_resolved$location_family
        .family_id_k_vector <- .het_resolved$family_id_k_vector
      } else {
        family_promoted     <- .gdpar_promote_scope_per_observation(
          family, names(amm_list_canonical)
        )
        .family_id_k_vector <- NULL
      }
      return(.gdpar_K(
        amm_list_canonical = amm_list_canonical,
        family             = family_promoted,
        data               = data,
        prior              = prior,
        anchor             = anchor,
        outcome_name       = outcome_name,
        formula_env        = formula_env,
        skip_id_check      = skip_id_check,
        chains             = chains,
        iter_warmup        = iter_warmup,
        iter_sampling      = iter_sampling,
        adapt_delta        = adapt_delta,
        max_treedepth      = max_treedepth,
        refresh            = refresh,
        verbose            = verbose,
        seed               = seed,
        group              = group,
        parametrization    = parametrization,
        parametrization_a  = parametrization_a,
        parametrization_W  = parametrization_W,
        id_check_rigor     = id_check_rigor,
        family_id_k_vector = .family_id_k_vector,
        call               = call,
        ...
      ))
    }
    if (.family_is_named_list) {
      gdpar_abort(
        paste(
          "A named list was supplied as 'family' (heterogeneous families per",
          "slot path, sub-phase 8.3.7) but the resolved K equals 1. The",
          "heterogeneous path requires K >= 2 K-individual slots. Pass a",
          "single gdpar_family object for the K = 1 regime."
        ),
        class = "gdpar_input_error",
        data = list(K = K)
      )
    }
    k_name <- names(amm_list_canonical)[[1L]]
    amm <- amm_list_canonical[[1L]]
    vars_for_rhs <- union(
      if (!is.null(amm$a)) all.vars(amm$a) else character(0L),
      if (!is.null(amm$b)) all.vars(amm$b) else character(0L)
    )
    rhs_str <- if (length(vars_for_rhs) == 0L) "1" else
               paste(vars_for_rhs, collapse = " + ")
    formula <- stats::as.formula(
      paste(outcome_name, "~", rhs_str),
      env = formula_env
    )
    family <- .gdpar_promote_scope_per_observation(family, k_name)
  } else if (!is.null(W)) {
    gdpar_abort(
      paste(
        "The 'W' argument of gdpar() is reserved for the K-individual paths",
        "(formula set, named list of amm_spec, or classic formula with AMM",
        "wrappers in the RHS). In the legacy single-amm_spec path, the",
        "modulating basis travels via amm$W (see amm_spec(W = W_basis(...)))."
      ),
      class = "gdpar_input_error"
    )
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    gdpar_abort(
      "Argument 'formula' must be a two-sided formula 'y ~ ...'.",
      class = "gdpar_input_error"
    )
  }
  assert_inherits(amm, "amm_spec", "amm")
  if (inherits(family, "gdpar_family_multi") &&
      (is.null(amm$p) || amm$p == 1L)) {
    gdpar_abort(
      paste(
        "A multivariate family (gdpar_family_multi) was supplied but",
        "amm$p is NULL or 1. Use amm_spec(p = ...) with p > 1 to match",
        "the family, or pass a univariate gdpar_family object."
      ),
      class = "gdpar_input_error",
      data = list(family_class = class(family), amm_p = amm$p)
    )
  }
  if (!is.null(amm$p) && amm$p > 1L) {
    return(.gdpar_multi(
      formula = formula, family = family, amm = amm, data = data,
      prior = prior, anchor = anchor, skip_id_check = skip_id_check,
      chains = chains, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, adapt_delta = adapt_delta,
      max_treedepth = max_treedepth, refresh = refresh,
      verbose = verbose, seed = seed,
      group = group,
      parametrization = parametrization,
      parametrization_a = parametrization_a,
      parametrization_W = parametrization_W,
      parametrization_aggregation = parametrization_aggregation,
      id_check_rigor = id_check_rigor,
      call = call, ...
    ))
  }
  assert_data_frame(data, "data")
  if (is.null(prior)) {
    prior <- gdpar_prior()
  }
  assert_inherits(prior, "gdpar_prior", "prior")
  assert_count(chains, "chains")
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_numeric_scalar(adapt_delta, "adapt_delta",
                        lower = 0.5, upper = 0.999)
  assert_count(max_treedepth, "max_treedepth")
  if (!is.numeric(refresh) || length(refresh) != 1L || refresh < 0) {
    gdpar_abort(
      "Argument 'refresh' must be a non-negative integer scalar.",
      class = "gdpar_input_error"
    )
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort(
      "Argument 'verbose' must be a logical scalar.",
      class = "gdpar_input_error"
    )
  }
  require_suggested("cmdstanr",
                    "fit Path 1 hierarchical Bayesian models")
  outcome_name <- as.character(formula[[2]])
  if (!outcome_name %in% colnames(data)) {
    gdpar_abort(
      sprintf(
        "Outcome variable '%s' not found in data.", outcome_name
      ),
      class = "gdpar_input_error"
    )
  }
  y <- data[[outcome_name]]
  y_bad <- if (is.numeric(y)) !is.finite(y) else is.na(y)
  if (any(y_bad)) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' contains %d non-finite value(s) (NA, NaN or Inf). Path 1 does not impute; please remove or fix them before fitting.",
        outcome_name, sum(y_bad)
      ),
      class = "gdpar_input_error"
    )
  }
  rhs <- formula[c(1L, 3L)]
  rhs <- stats::update(rhs, ~ . + 0)
  if (!is.null(amm$W)) {
    amm$W <- materialize_W_basis(amm$W, p = 1L)
  }
  design <- build_amm_design(amm, data, formula_rhs = rhs)
  anchor_value <- resolve_anchor(anchor, family, y, prior, verbose)
  id_report <- if (skip_id_check) {
    if (verbose) {
      gdpar_inform(
        "Skipping basis-restricted identifiability check (skip_id_check = TRUE).",
        class = "gdpar_identifiability_message"
      )
    }
    NULL
  } else {
    diag_theta <- if (!is.null(amm$b) && abs(anchor_value) < 1e-8) {
      1
    } else {
      anchor_value
    }
    rep <- gdpar_check_identifiability(
      amm = amm, data = data,
      theta_ref_init = diag_theta,
      formula_rhs = rhs
    )
    if (!isTRUE(rep$passed)) {
      gdpar_abort(
        paste(
          "Basis-restricted identifiability check failed at the diagnostic",
          "test point. Inspect the returned report or pass",
          "skip_id_check = TRUE to override at your own risk."
        ),
        class = "gdpar_identifiability_error",
        data = list(report = rep)
      )
    }
    rep
  }
  if (family$did_status == "holds_under_condition" && verbose) {
    gdpar_inform(
      sprintf(
        paste(
          "Family '%s' has D-ID status 'holds_under_condition'.",
          "Condition: %s. Reference: %s. The package documents this",
          "but does not verify it from data."
        ),
        family$name, family$did_condition, family$did_reference
      ),
      class = "gdpar_did_message"
    )
  }
  group_info <- .resolve_group_argument(group, data, n = length(y),
                                        verbose = verbose)
  if (!is.null(group_info)) {
    .check_group_aliasing_c7(
      design         = design,
      group_id       = group_info$group_id,
      group_var_name = group_info$var_name
    )
  }
  stan_data <- assemble_stan_data(
    design, family, amm, y, anchor_value,
    group_id = if (is.null(group_info)) NULL else group_info$group_id
  )
  parametrization_resolved <- resolve_parametrization(
    parametrization      = parametrization,
    parametrization_a    = parametrization_a,
    parametrization_W    = parametrization_W,
    prior                = prior,
    stan_data            = stan_data,
    amm                  = amm,
    preflight_seed       = seed,
    verbose              = verbose
  )
  stan_src <- generate_stan_code(
    prior,
    cp_a = parametrization_resolved$cp_a,
    cp_W = parametrization_resolved$cp_W
  )
  stan_path <- write_stan_to_tempfile(stan_src)
  cs_model <- cmdstanr::cmdstan_model(stan_path)
  sample_args <- list(
    data            = stan_data,
    chains          = as.integer(chains),
    iter_warmup     = as.integer(iter_warmup),
    iter_sampling   = as.integer(iter_sampling),
    adapt_delta     = adapt_delta,
    max_treedepth   = as.integer(max_treedepth),
    refresh         = as.integer(refresh),
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) {
    sample_args$seed <- as.integer(seed)
  }
  extra <- list(...)
  for (nm in names(extra)) sample_args[[nm]] <- extra[[nm]]
  fit <- if (isTRUE(verbose)) {
    do.call(cs_model$sample, sample_args)
  } else {
    msg_file <- tempfile(fileext = ".log")
    msg_con <- file(msg_file, open = "wt")
    sink(msg_con, type = "message")
    flush_captured <- function() {
      if (sink.number(type = "message") > 0L) {
        sink(NULL, type = "message")
      }
      if (isOpen(msg_con)) close(msg_con)
      captured <- tryCatch(readLines(msg_file, warn = FALSE),
                           error = function(e) character(0))
      unlink(msg_file)
      deduped <- dedup_message_blocks(captured)
      for (line in deduped) message(line)
    }
    on.exit(flush_captured(), add = TRUE, after = FALSE)
    do.call(cs_model$sample, sample_args)
  }
  diagnostics <- compute_diagnostics(fit, verbose = verbose)
  obj <- list(
    fit = fit,
    amm = amm,
    family = family,
    prior = prior,
    design = design,
    anchor = anchor_value,
    stan_data = stan_data,
    identifiability_report = id_report,
    diagnostics = diagnostics,
    parametrization = list(
      cp_a = parametrization_resolved$cp_a,
      cp_W = parametrization_resolved$cp_W,
      meta = parametrization_resolved$meta
    ),
    group_info = group_info,
    call = call,
    path = path
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

#' Internal: multivariate path of gdpar() for amm$p > 1
#'
#' Mirror of the scalar path, with three structural differences: (a)
#' the outcome is a matrix of dimensions n by p stored as a single
#' matrix-valued column of \code{data}; (b) the family is automatically
#' promoted to \code{gdpar_family_multi} when the user supplies a
#' univariate \code{gdpar_family}; (c) Stan codegen, design assembly,
#' and assembler dispatch to the multivariate helpers added in Phase F.
#' Pre-flight CP/NCP per-coordinate is deferred to Phase H of Block
#' 5.2; the parametrization arguments are honored literally and the
#' helper rejects \code{cp_a = TRUE} with an explicit deferred-feature
#' error.
#'
#' @keywords internal
#' @noRd
.gdpar_multi <- function(formula, family, amm, data, prior, anchor,
                         skip_id_check, chains, iter_warmup,
                         iter_sampling, adapt_delta, max_treedepth,
                         refresh, verbose, seed,
                         group = NULL,
                         parametrization,
                         parametrization_a, parametrization_W,
                         parametrization_aggregation = NULL,
                         id_check_rigor = "full",
                         call, ...) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    gdpar_abort(
      "Argument 'formula' must be a two-sided formula 'y ~ ...'.",
      class = "gdpar_input_error"
    )
  }
  assert_data_frame(data, "data")
  if (is.null(prior)) {
    prior <- gdpar_prior()
  }
  assert_inherits(prior, "gdpar_prior", "prior")
  assert_count(chains, "chains")
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_numeric_scalar(adapt_delta, "adapt_delta",
                        lower = 0.5, upper = 0.999)
  assert_count(max_treedepth, "max_treedepth")
  if (!is.numeric(refresh) || length(refresh) != 1L || refresh < 0) {
    gdpar_abort(
      "Argument 'refresh' must be a non-negative integer scalar.",
      class = "gdpar_input_error"
    )
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort(
      "Argument 'verbose' must be a logical scalar.",
      class = "gdpar_input_error"
    )
  }
  require_suggested("cmdstanr",
                    "fit Path 1 hierarchical Bayesian models")
  p <- amm$p
  if (inherits(family, "gdpar_family") &&
      !inherits(family, "gdpar_family_multi")) {
    if (isTRUE(verbose)) {
      gdpar_inform(
        sprintf(
          paste(
            "Auto-promoted univariate family '%s' to gdpar_family_multi",
            "via gdpar_family_multi(family, p = %d) because amm$p > 1.",
            "Pass an explicit gdpar_family_multi to silence."
          ),
          family$name, p
        ),
        class = "gdpar_family_promotion_message"
      )
    }
    family <- gdpar_family_multi(family, p = p)
  }
  if (!inherits(family, "gdpar_family_multi")) {
    gdpar_abort(
      paste(
        "When amm$p > 1L, 'family' must be of class gdpar_family",
        "(auto-promoted) or gdpar_family_multi. Received:",
        paste(class(family), collapse = ", "), "."
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(family))
    )
  }
  if (family$p != p) {
    gdpar_abort(
      sprintf(
        "Family p = %d does not match amm$p = %d. They must coincide.",
        family$p, p
      ),
      class = "gdpar_input_error",
      data = list(family_p = family$p, amm_p = p)
    )
  }
  outcome_name <- as.character(formula[[2]])
  if (!outcome_name %in% colnames(data)) {
    gdpar_abort(
      sprintf("Outcome variable '%s' not found in data.", outcome_name),
      class = "gdpar_input_error"
    )
  }
  y <- data[[outcome_name]]
  if (!is.matrix(y) && !is.array(y)) {
    gdpar_abort(
      sprintf(
        paste(
          "When amm$p > 1L (multivariate path), the outcome '%s' in",
          "'data' must be a matrix of dimensions n by p = %d. Store the",
          "outcome as a matrix column with data$%s <- cbind(...) or as",
          "a numeric matrix and pass it via the data frame."
        ),
        outcome_name, p, outcome_name
      ),
      class = "gdpar_input_error",
      data = list(outcome_class = class(y))
    )
  }
  y_bad <- if (is.numeric(y)) !is.finite(y) else is.na(y)
  if (any(y_bad)) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' contains %d non-finite value(s) (NA, NaN or Inf). Path 1 does not impute; please remove or fix them before fitting.",
        outcome_name, sum(y_bad)
      ),
      class = "gdpar_input_error"
    )
  }
  if (ncol(y) != p) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' has %d columns but amm$p = %d. They must coincide.",
        outcome_name, ncol(y), p
      ),
      class = "gdpar_input_error",
      data = list(outcome_cols = ncol(y), p = p)
    )
  }
  rhs <- formula[c(1L, 3L)]
  rhs <- stats::update(rhs, ~ . + 0)
  if (!is.null(amm$W)) {
    amm$W <- materialize_W_basis(amm$W, p = p)
  }
  design <- build_amm_design(amm, data, formula_rhs = rhs)
  anchor_value <- resolve_anchor_multi(anchor, family, y, prior, p,
                                       verbose)
  id_report <- if (skip_id_check) {
    if (verbose) {
      gdpar_inform(
        "Skipping basis-restricted identifiability check (skip_id_check = TRUE).",
        class = "gdpar_identifiability_message"
      )
    }
    NULL
  } else {
    diag_theta_multi <- if (!is.null(amm$b) && all(abs(anchor_value) < 1e-8)) {
      rep(1, p)
    } else {
      anchor_value
    }
    rep_multi <- gdpar_check_identifiability(
      amm            = amm,
      data           = data,
      theta_ref_init = diag_theta_multi,
      formula_rhs    = rhs,
      rigor          = id_check_rigor
    )
    if (!isTRUE(rep_multi$passed)) {
      gdpar_abort(
        paste(
          "Basis-restricted identifiability check (C1-C4 + C4-bis)",
          "failed for the multivariate spec at the diagnostic test",
          "point. Inspect the returned report or pass skip_id_check =",
          "TRUE to override at your own risk."
        ),
        class = "gdpar_identifiability_error",
        data = list(report = rep_multi)
      )
    }
    rep_multi
  }
  if (family$did_status == "holds_under_condition" && verbose) {
    gdpar_inform(
      sprintf(
        paste(
          "Family '%s' has D-ID status 'holds_under_condition'.",
          "Condition: %s. Reference: %s. The package documents this",
          "but does not verify it from data."
        ),
        family$name, family$did_condition, family$did_reference
      ),
      class = "gdpar_did_message"
    )
  }
  group_info <- .resolve_group_argument(group, data, n = nrow(y),
                                        verbose = verbose)
  if (!is.null(group_info)) {
    .check_group_aliasing_c7(
      design         = design,
      group_id       = group_info$group_id,
      group_var_name = group_info$var_name
    )
  }
  stan_data <- assemble_stan_data(
    design, family, amm, y, anchor_value,
    group_id = if (is.null(group_info)) NULL else group_info$group_id
  )
  resolved_multi <- resolve_parametrization_multi(
    parametrization             = parametrization,
    parametrization_a           = parametrization_a,
    parametrization_W           = parametrization_W,
    parametrization_aggregation = parametrization_aggregation,
    prior                       = prior,
    stan_data                   = stan_data,
    amm                         = amm,
    preflight_seed              = seed,
    verbose                     = verbose
  )
  parametrization_resolved <- list(
    cp_a         = resolved_multi$cp_a,
    cp_W         = resolved_multi$cp_W,
    cp_a_per_k   = resolved_multi$cp_a_per_k,
    cp_W_per_k   = resolved_multi$cp_W_per_k,
    report       = resolved_multi$report,
    meta = list(
      mode = "multivariate_phase_h",
      note = paste(
        "Phase H.2: preflight CP/NCP per-k (Path B' per coord) and",
        "aggregation via parametrization_aggregation (default",
        "'any_ncp'); per-k Stan wiring active via segment()-based",
        "per-coordinate priors. Uniform-NCP and uniform-CP branches",
        "compile to byte-identical code w.r.t. the H.1 template."
      ),
      aggregation = resolved_multi$aggregation,
      requested = list(
        parametrization             = parametrization,
        parametrization_a           = parametrization_a,
        parametrization_W           = parametrization_W,
        parametrization_aggregation = parametrization_aggregation
      )
    )
  )
  stan_src <- generate_stan_code_multi(
    prior,
    cp_a       = isTRUE(parametrization_resolved$cp_a),
    cp_W       = isTRUE(parametrization_resolved$cp_W),
    cp_a_per_k = parametrization_resolved$cp_a_per_k
  )
  if (length(unique(parametrization_resolved$cp_a_per_k)) > 1L) {
    stan_data$cp_a_per_k_data <- as.integer(
      parametrization_resolved$cp_a_per_k
    )
  }
  stan_path <- write_stan_to_tempfile(stan_src)
  cs_model <- cmdstanr::cmdstan_model(stan_path)
  sample_args <- list(
    data            = stan_data,
    chains          = as.integer(chains),
    iter_warmup     = as.integer(iter_warmup),
    iter_sampling   = as.integer(iter_sampling),
    adapt_delta     = adapt_delta,
    max_treedepth   = as.integer(max_treedepth),
    refresh         = as.integer(refresh),
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) {
    sample_args$seed <- as.integer(seed)
  }
  extra <- list(...)
  for (nm in names(extra)) sample_args[[nm]] <- extra[[nm]]
  fit <- if (isTRUE(verbose)) {
    do.call(cs_model$sample, sample_args)
  } else {
    msg_file <- tempfile(fileext = ".log")
    msg_con <- file(msg_file, open = "wt")
    sink(msg_con, type = "message")
    flush_captured <- function() {
      if (sink.number(type = "message") > 0L) {
        sink(NULL, type = "message")
      }
      if (isOpen(msg_con)) close(msg_con)
      captured <- tryCatch(readLines(msg_file, warn = FALSE),
                           error = function(e) character(0))
      unlink(msg_file)
      deduped <- dedup_message_blocks(captured)
      for (line in deduped) message(line)
    }
    on.exit(flush_captured(), add = TRUE, after = FALSE)
    do.call(cs_model$sample, sample_args)
  }
  diagnostics <- compute_diagnostics(fit, verbose = verbose)
  obj <- list(
    fit = fit,
    amm = amm,
    family = family,
    prior = prior,
    design = design,
    anchor = anchor_value,
    stan_data = stan_data,
    identifiability_report = id_report,
    diagnostics = diagnostics,
    parametrization = parametrization_resolved,
    group_info = group_info,
    call = call,
    path = "bayes",
    p = p
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

#' Internal: K-individual path of gdpar() for K > 1 with p = 1
#'
#' Mirror of the K = 1 path tailored to the distributional regression
#' template \code{amm_distrib_K.stan} (Sub-phase 8.3.3 Unit 3,
#' 2026-05-20). Three structural differences relative to
#' \code{.gdpar_multi()}:
#' \enumerate{
#'   \item The outcome \code{y} is a length-n vector (the K slots all
#'     parameterize \emph{one} univariate response via canonical
#'     links).
#'   \item The family is a univariate \code{gdpar_family} whose
#'     K-individual slots have been promoted to
#'     \code{scope = "per_observation"} via
#'     \code{.gdpar_promote_scope_per_observation()} by the caller.
#'   \item Stan codegen, design assembly and the data assembler
#'     dispatch to \code{.build_amm_design_K()},
#'     \code{.assemble_stan_data_K()} and \code{generate_stan_code_K()}
#'     respectively, all of which target the new template
#'     \code{amm_distrib_K.stan}.
#' }
#'
#' Pre-flight CP/NCP per-slot is not wired in Unit 3 (deuda potencial
#' for sub-phase 8.3.4 onwards): the helper honors \code{parametrization}
#' as a uniform flag across slots and falls back to NCP when
#' \code{parametrization = "auto"}.
#'
#' @keywords internal
#' @noRd
.gdpar_K <- function(amm_list_canonical, family, data, prior, anchor,
                     outcome_name, formula_env,
                     skip_id_check, chains, iter_warmup,
                     iter_sampling, adapt_delta, max_treedepth,
                     refresh, verbose, seed,
                     group = NULL,
                     parametrization,
                     parametrization_a, parametrization_W,
                     id_check_rigor = "full",
                     family_id_k_vector = NULL,
                     call, ...) {
  assert_data_frame(data, "data")
  if (is.null(prior)) {
    prior <- gdpar_prior()
  }
  assert_inherits(prior, "gdpar_prior", "prior")
  assert_count(chains, "chains")
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_numeric_scalar(adapt_delta, "adapt_delta",
                        lower = 0.5, upper = 0.999)
  assert_count(max_treedepth, "max_treedepth")
  if (!is.numeric(refresh) || length(refresh) != 1L || refresh < 0) {
    gdpar_abort(
      "Argument 'refresh' must be a non-negative integer scalar.",
      class = "gdpar_input_error"
    )
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort(
      "Argument 'verbose' must be a logical scalar.",
      class = "gdpar_input_error"
    )
  }
  require_suggested("cmdstanr",
                    "fit Path 1 hierarchical Bayesian models")
  bd <- .gdpar_K_build(
    amm_list_canonical = amm_list_canonical, family = family,
    data = data, prior = prior, anchor = anchor,
    outcome_name = outcome_name, formula_env = formula_env,
    skip_id_check = skip_id_check, verbose = verbose, group = group,
    parametrization = parametrization,
    parametrization_a = parametrization_a,
    parametrization_W = parametrization_W,
    id_check_rigor = id_check_rigor,
    family_id_k_vector = family_id_k_vector,
    compile_model_methods = FALSE
  )
  cs_model <- bd$cs_model
  stan_data <- bd$stan_data
  amm_list_canonical <- bd$meta$amm_list_canonical
  family <- bd$meta$family
  prior <- bd$meta$prior
  design_K <- bd$meta$design_K
  anchor_value <- bd$meta$anchor_value
  id_report <- bd$meta$id_report
  group_info <- bd$meta$group_info
  parametrization_resolved <- bd$meta$parametrization_resolved
  slot_names <- bd$meta$slot_names
  K <- bd$meta$K
  sample_args <- list(
    data            = stan_data,
    chains          = as.integer(chains),
    iter_warmup     = as.integer(iter_warmup),
    iter_sampling   = as.integer(iter_sampling),
    adapt_delta     = adapt_delta,
    max_treedepth   = as.integer(max_treedepth),
    refresh         = as.integer(refresh),
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) {
    sample_args$seed <- as.integer(seed)
  }
  if (!is.null(bd$init)) {
    sample_args$init <- bd$init
  }
  extra <- list(...)
  for (nm in names(extra)) sample_args[[nm]] <- extra[[nm]]
  fit <- if (isTRUE(verbose)) {
    do.call(cs_model$sample, sample_args)
  } else {
    msg_file <- tempfile(fileext = ".log")
    msg_con <- file(msg_file, open = "wt")
    sink(msg_con, type = "message")
    flush_captured <- function() {
      if (sink.number(type = "message") > 0L) {
        sink(NULL, type = "message")
      }
      if (isOpen(msg_con)) close(msg_con)
      captured <- tryCatch(readLines(msg_file, warn = FALSE),
                           error = function(e) character(0))
      unlink(msg_file)
      deduped <- dedup_message_blocks(captured)
      for (line in deduped) message(line)
    }
    on.exit(flush_captured(), add = TRUE, after = FALSE)
    do.call(cs_model$sample, sample_args)
  }
  diagnostics <- compute_diagnostics(fit, verbose = verbose)
  use_groups_int <- if (is.null(group_info)) 0L else 1L
  did_post_fit <- .compute_info_ratio_K(
    fit         = fit,
    family      = family,
    slot_names  = slot_names,
    use_groups  = use_groups_int,
    prior       = prior
  )
  if (isTRUE(did_post_fit$any_info_error)) {
    info_err_slots <- vapply(did_post_fit$per_slot,
      function(s) if (identical(s$status, "information_error")) s$slot
                  else NA_character_,
      character(1L))
    info_err_slots <- info_err_slots[!is.na(info_err_slots)]
    gdpar_warn(
      sprintf(
        paste(
          "Post-fit information contraction below threshold 0.1 for",
          "slot(s) %s: the posterior recovers the prior almost",
          "exactly, indicating that the data is essentially",
          "uninformative for these slot anchor(s). Inspect",
          "$identifiability_post_fit for per-slot diagnostics."
        ),
        paste(sQuote(info_err_slots), collapse = ", ")
      ),
      class = "gdpar_information_error",
      data  = list(report = did_post_fit, slots = info_err_slots)
    )
  } else if (isTRUE(did_post_fit$any_warn) && verbose) {
    warn_slots <- vapply(did_post_fit$per_slot,
      function(s) if (identical(s$status, "warn")) s$slot
                  else NA_character_,
      character(1L))
    warn_slots <- warn_slots[!is.na(warn_slots)]
    gdpar_warn(
      sprintf(
        paste(
          "Post-fit information contraction in [0.1, 0.5) for slot(s)",
          "%s: weak learning. Inspect $identifiability_post_fit."
        ),
        paste(sQuote(warn_slots), collapse = ", ")
      ),
      class = "gdpar_information_warning",
      data  = list(report = did_post_fit, slots = warn_slots)
    )
  }
  obj <- list(
    fit                       = fit,
    amm_list_canonical        = amm_list_canonical,
    family                    = family,
    prior                     = prior,
    design_K                  = design_K,
    anchor                    = anchor_value,
    stan_data                 = stan_data,
    identifiability_report    = id_report,
    identifiability_post_fit  = did_post_fit,
    diagnostics               = diagnostics,
    parametrization           = parametrization_resolved,
    group_info                = group_info,
    call                      = call,
    path                      = "bayes",
    K                         = K,
    slot_names                = slot_names
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

#' Internal: assemble the K-individual cmdstan model and data (pre-sampling)
#'
#' Behaviour-preserving extraction of the build phase of
#' \code{.gdpar_K()} (everything up to but not including the call to
#' \code{cs_model$sample}). \code{.gdpar_K()} and the opt-in
#' geometry-adaptive entry point \code{gdpar_geom_fit()} share this
#' single source so the model + data are built and compiled exactly once
#' (no duplication, no throwaway computation). The default branch
#' (\code{compile_model_methods = FALSE}) compiles the model with the
#' identical call \code{cmdstanr::cmdstan_model(stan_path)} that
#' \code{.gdpar_K()} used before the refactor, so the K-path goldens
#' remain bit-identical; \code{gdpar_geom_fit()} passes
#' \code{compile_model_methods = TRUE} to expose
#' \code{$log_prob}/\code{$grad_log_prob}/\code{$hessian} for the
#' geometry engine. Canonised in Block RG, RG.6 part ii (session B9.31).
#'
#' @return A list with \code{cs_model} (the compiled
#'   \code{CmdStanModel}), \code{stan_data}, \code{init} (the per-chain
#'   init function for the Tweedie \code{p} slot, or \code{NULL}) and
#'   \code{meta} (every object \code{.gdpar_K()} needs to assemble the
#'   final \code{gdpar_fit} after sampling).
#' @keywords internal
#' @noRd
.gdpar_K_build <- function(amm_list_canonical, family, data, prior, anchor,
                           outcome_name, formula_env, skip_id_check, verbose,
                           group = NULL, parametrization,
                           parametrization_a, parametrization_W,
                           id_check_rigor = "full",
                           family_id_k_vector = NULL,
                           compile_model_methods = FALSE) {
  if (is.null(prior)) {
    prior <- gdpar_prior()
  }
  require_suggested("cmdstanr",
                    "fit Path 1 hierarchical Bayesian models")
  slot_names <- names(amm_list_canonical)
  K <- length(amm_list_canonical)
  if (!outcome_name %in% colnames(data)) {
    gdpar_abort(
      sprintf("Outcome variable '%s' not found in data.", outcome_name),
      class = "gdpar_input_error"
    )
  }
  y <- data[[outcome_name]]
  if (is.matrix(y) || is.array(y) && length(dim(y)) > 1L) {
    gdpar_abort(
      sprintf(
        paste(
          "When K > 1 with p = 1, outcome '%s' must be a length-n vector",
          "(one univariate response shared across the K slots). For",
          "matrix-valued outcomes use the K = 1 + p > 1 multi-coord path."
        ),
        outcome_name
      ),
      class = "gdpar_input_error",
      data = list(outcome_class = class(y))
    )
  }
  y_bad <- if (is.numeric(y)) !is.finite(y) else is.na(y)
  if (any(y_bad)) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' contains %d non-finite value(s) (NA, NaN or Inf). Path 1 does not impute; please remove or fix them before fitting.",
        outcome_name, sum(y_bad)
      ),
      class = "gdpar_input_error"
    )
  }
  union_vars <- character(0)
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]
    if (!is.null(a_k$a)) union_vars <- c(union_vars, all.vars(a_k$a))
    if (!is.null(a_k$b)) union_vars <- c(union_vars, all.vars(a_k$b))
  }
  union_vars <- unique(union_vars)
  rhs_str <- if (length(union_vars) == 0L) "1" else
             paste(union_vars, collapse = " + ")
  formula <- stats::as.formula(
    paste(outcome_name, "~", rhs_str),
    env = formula_env
  )
  rhs <- formula[c(1L, 3L)]
  rhs <- stats::update(rhs, ~ . + 0)
  for (k in seq_len(K)) {
    if (!is.null(amm_list_canonical[[k]]$W)) {
      amm_list_canonical[[k]]$W <- materialize_W_basis(
        amm_list_canonical[[k]]$W, p = 1L
      )
    }
  }
  design_K <- .build_amm_design_K(amm_list_canonical, data,
                                   formula_rhs = rhs)
  anchor_value <- resolve_anchor_K(anchor, family, y, prior,
                                   slot_names, verbose)
  id_report <- if (skip_id_check) {
    if (verbose) {
      gdpar_inform(
        "Skipping basis-restricted identifiability check (skip_id_check = TRUE).",
        class = "gdpar_identifiability_message"
      )
    }
    NULL
  } else {
    per_slot_reports <- vector("list", K)
    names(per_slot_reports) <- slot_names
    for (k in seq_len(K)) {
      a_k <- amm_list_canonical[[k]]
      diag_theta <- if (!is.null(a_k$b) && abs(anchor_value[[k]]) < 1e-8) {
        1
      } else {
        anchor_value[[k]]
      }
      slot_vars <- unique(c(
        if (!is.null(a_k$a)) all.vars(a_k$a) else character(0),
        if (!is.null(a_k$b)) all.vars(a_k$b) else character(0)
      ))
      slot_rhs_str <- if (length(slot_vars) == 0L) "1" else
                       paste(slot_vars, collapse = " + ")
      slot_rhs <- stats::as.formula(
        paste("~", slot_rhs_str), env = formula_env
      )
      rep_k <- gdpar_check_identifiability(
        amm            = a_k, data = data,
        theta_ref_init = diag_theta,
        formula_rhs    = slot_rhs
      )
      if (!isTRUE(rep_k$passed)) {
        gdpar_abort(
          sprintf(
            paste(
              "Basis-restricted identifiability check failed for slot '%s' at",
              "the diagnostic test point. Inspect the returned report or pass",
              "skip_id_check = TRUE to override at your own risk."
            ),
            slot_names[k]
          ),
          class = "gdpar_identifiability_error",
          data = list(slot = slot_names[k], report = rep_k)
        )
      }
      per_slot_reports[[k]] <- rep_k
    }
    K_level <- .check_identifiability_K(design_K,
                                         rigor = id_check_rigor)
    if (!isTRUE(K_level$passed)) {
      failed_part <- if (!K_level$per_slot_rank$passed) {
        "per-slot Z_a^{(k)} rank check (D-B3)"
      } else {
        "cross-slot extended Gram check (D-B2)"
      }
      gdpar_abort(
        sprintf(
          paste(
            "K-level identifiability check failed at the %s. Inspect",
            "the returned report or pass skip_id_check = TRUE to",
            "override at your own risk."
          ),
          failed_part
        ),
        class = "gdpar_identifiability_error",
        data = list(report = K_level)
      )
    }
    attr(per_slot_reports, "K_level") <- K_level
    per_slot_reports
  }
  if (family$did_status == "holds_under_condition" && verbose) {
    gdpar_inform(
      sprintf(
        paste(
          "Family '%s' has D-ID status 'holds_under_condition'.",
          "Condition: %s. Reference: %s. The package documents this",
          "but does not verify it from data."
        ),
        family$name, family$did_condition, family$did_reference
      ),
      class = "gdpar_did_message"
    )
  }
  group_info <- .resolve_group_argument(group, data, n = length(y),
                                        verbose = verbose)
  if (!is.null(group_info)) {
    for (k in seq_len(K)) {
      slot_design <- list(
        Z_a = design_K$Z_a_k_list[[k]],
        Z_b = design_K$Z_b_k_list[[k]],
        X   = design_K$X
      )
      .check_group_aliasing_c7(
        design         = slot_design,
        group_id       = group_info$group_id,
        group_var_name = group_info$var_name
      )
    }
  }
  stan_data <- .assemble_stan_data_K(
    design_K, family, amm_list_canonical, y, anchor_value,
    group_id = if (is.null(group_info)) NULL else group_info$group_id,
    family_id_k_vector = family_id_k_vector
  )
  cp_a <- if (!is.null(parametrization_a)) {
    identical(parametrization_a, "cp")
  } else {
    identical(parametrization, "cp")
  }
  cp_W <- if (!is.null(parametrization_W)) {
    identical(parametrization_W, "cp")
  } else {
    identical(parametrization, "cp")
  }
  parametrization_resolved <- list(
    cp_a       = isTRUE(cp_a),
    cp_W       = isTRUE(cp_W),
    cp_a_per_K = NULL,
    meta = list(
      mode = "K_individual_unit3",
      note = paste(
        "Unit 3 (8.3.3) ships with uniform CP/NCP across the K slots;",
        "per-slot preflight (cp_a_per_K) is queued for sub-phase 8.3.4",
        "or later."
      ),
      requested = list(
        parametrization   = parametrization,
        parametrization_a = parametrization_a,
        parametrization_W = parametrization_W
      )
    )
  )
  stan_src <- generate_stan_code_K(
    prior,
    cp_a       = parametrization_resolved$cp_a,
    cp_W       = parametrization_resolved$cp_W,
    cp_a_per_K = parametrization_resolved$cp_a_per_K,
    family     = family
  )
  stan_path <- write_stan_to_tempfile(stan_src)
  cs_model <- if (isTRUE(compile_model_methods)) {
    cmdstanr::cmdstan_model(stan_path, compile_model_methods = TRUE)
  } else {
    cmdstanr::cmdstan_model(stan_path)
  }
  # Tweedie p slot carries the structural constraint p in (1.01, 1.99)
  # via the uniform prior. The default cmdstanr init draws unconstrained
  # parameters in (-2, 2), which lands outside the support often enough
  # to abort all chains with -inf log_prob; pin theta_ref_k slot p at the
  # midpoint 1.5 per chain. Identical to the previous inline .gdpar_K init.
  init_fn <- NULL
  if (!is.null(family$stan_id) && !is.na(family$stan_id) &&
      as.integer(family$stan_id) == 9L) {
    K_loc <- as.integer(stan_data$K)
    J_groups_loc <- max(1L, as.integer(stan_data$J_groups))
    init_theta_ref <- matrix(0, nrow = J_groups_loc, ncol = K_loc)
    init_theta_ref[, K_loc] <- 1.5
    init_fn <- function() list(theta_ref_k = init_theta_ref)
  }
  list(
    cs_model  = cs_model,
    stan_data = stan_data,
    init      = init_fn,
    meta = list(
      amm_list_canonical       = amm_list_canonical,
      family                   = family,
      prior                    = prior,
      design_K                 = design_K,
      anchor_value             = anchor_value,
      id_report                = id_report,
      group_info               = group_info,
      parametrization_resolved = parametrization_resolved,
      slot_names               = slot_names,
      K                        = K,
      stan_src                 = stan_src
    )
  )
}

#' Internal: resolve the anchor argument for the K-individual path
#'
#' K-individual counterpart of \code{resolve_anchor} and
#' \code{resolve_anchor_multi}. Accepts a numeric scalar (broadcast to
#' length K), a numeric vector of length K (optionally named to enforce
#' alignment with \code{slot_names}), \code{"prior_mean"} (vector of
#' zeros), or \code{"empirical_y"} (location slot anchored at the
#' link-transformed mean of \code{y}, every other slot anchored at
#' zero; the per-slot link semantics for dispersion slots are deferred
#' to sub-phase 8.3.4 when Beta and Gamma K = 2 introduce richer link
#' choices).
#'
#' @param anchor Either a numeric scalar (broadcast), a numeric vector
#'   of length K (optionally named), or one of \code{"prior_mean"},
#'   \code{"empirical_y"}.
#' @param family A \code{gdpar_family} object whose K-individual slots
#'   have been promoted to per-observation scope.
#' @param y Numeric vector of outcomes (length n).
#' @param prior A \code{gdpar_prior} object.
#' @param slot_names Character vector of length K with the canonical
#'   slot names of the K-individual specification.
#' @param verbose Logical; whether to inform about the chosen anchor.
#' @return Numeric vector of length K with names matching
#'   \code{slot_names}.
#' @keywords internal
#' @noRd
resolve_anchor_K <- function(anchor, family, y, prior, slot_names,
                             verbose) {
  K <- length(slot_names)
  if (is.numeric(anchor) && length(anchor) == 1L && is.finite(anchor)) {
    out <- rep(as.double(anchor), K)
    names(out) <- slot_names
    return(out)
  }
  if (is.numeric(anchor) && length(anchor) == K &&
      all(is.finite(anchor))) {
    if (!is.null(names(anchor)) && any(nzchar(names(anchor)))) {
      if (!setequal(names(anchor), slot_names)) {
        gdpar_abort(
          sprintf(
            paste(
              "Named anchor vector has names %s but the K-individual",
              "specification declares slots %s. These must coincide as sets."
            ),
            paste(sQuote(names(anchor)), collapse = ", "),
            paste(sQuote(slot_names), collapse = ", ")
          ),
          class = "gdpar_input_error",
          data = list(received_names = names(anchor),
                      expected_names = slot_names)
        )
      }
      out <- as.double(anchor[slot_names])
    } else {
      out <- as.double(anchor)
    }
    names(out) <- slot_names
    return(out)
  }
  if (is.character(anchor) && length(anchor) == 1L) {
    if (anchor == "prior_mean") {
      out <- rep(0.0, K)
      names(out) <- slot_names
      return(out)
    }
    if (anchor == "empirical_y") {
      yb <- mean(as.double(y))
      eta_loc <- tryCatch(
        family$linkfun(yb),
        error = function(e) {
          gdpar_abort(
            sprintf(
              "Failed to apply the location-slot link to the outcome mean: %s",
              conditionMessage(e)
            ),
            class = "gdpar_input_error"
          )
        }
      )
      out <- rep(0.0, K)
      out[1L] <- eta_loc
      names(out) <- slot_names
      if (verbose) {
        gdpar_inform(
          sprintf(
            paste(
              "Using empirical anchor on the location slot '%s' = link(mean(y))",
              "= %g; every other slot anchored at 0. Per-slot link semantics",
              "for dispersion slots land in sub-phase 8.3.4."
            ),
            slot_names[1L], eta_loc
          ),
          class = "gdpar_anchor_message"
        )
      }
      return(out)
    }
  }
  gdpar_abort(
    paste(
      "Argument 'anchor' must be a numeric scalar (broadcast to length K),",
      "a numeric vector of length K (optionally named with the slot names),",
      "'prior_mean', or 'empirical_y'."
    ),
    class = "gdpar_input_error",
    data = list(received = anchor, K = K, slot_names = slot_names)
  )
}

#' Internal: resolve the anchor argument for the multivariate path
#'
#' Multivariate counterpart of \code{resolve_anchor}. Accepts a numeric
#' scalar (broadcast to a vector of length \code{p}), a numeric vector
#' of length \code{p}, or one of the strings \code{"prior_mean"}
#' (vector of zeros), \code{"empirical_y"} (per-coordinate link applied
#' to the column means of the outcome matrix).
#'
#' @param anchor Either a numeric scalar, a numeric vector of length
#'   \code{p}, or one of \code{"prior_mean"}, \code{"empirical_y"}.
#' @param family A \code{gdpar_family_multi} object.
#' @param y Numeric matrix of outcomes (n by p).
#' @param prior A \code{gdpar_prior} object.
#' @param p Length of \code{theta_ref}.
#' @param verbose Logical; whether to inform about the chosen anchor.
#' @return Numeric vector of length \code{p}.
#' @keywords internal
#' @noRd
resolve_anchor_multi <- function(anchor, family, y, prior, p, verbose) {
  if (is.numeric(anchor) && length(anchor) == 1L && is.finite(anchor)) {
    return(rep(as.double(anchor), p))
  }
  if (is.numeric(anchor) && length(anchor) == p &&
      all(is.finite(anchor))) {
    return(as.double(anchor))
  }
  if (is.character(anchor) && length(anchor) == 1L) {
    if (anchor == "prior_mean") {
      return(rep(0.0, p))
    }
    if (anchor == "empirical_y") {
      yb <- colMeans(as.matrix(y))
      eta_anchor <- vapply(seq_len(p), function(k) {
        fam_k <- family$families[[k]]
        tryCatch(
          fam_k$linkfun(yb[k]),
          error = function(e) {
            gdpar_abort(
              sprintf(
                paste(
                  "Failed to apply the link function of coordinate %d",
                  "to the column mean of the outcome: %s"
                ),
                k, conditionMessage(e)
              ),
              class = "gdpar_input_error"
            )
          }
        )
      }, numeric(1L))
      if (verbose) {
        gdpar_inform(
          sprintf(
            "Using empirical anchor per coordinate: link(colMeans(y)) = [%s].",
            paste(format(eta_anchor, digits = 4), collapse = ", ")
          ),
          class = "gdpar_anchor_message"
        )
      }
      return(eta_anchor)
    }
  }
  gdpar_abort(
    paste(
      "Argument 'anchor' must be a numeric scalar (broadcast to length",
      "p), a numeric vector of length p, 'prior_mean', or",
      "'empirical_y'."
    ),
    class = "gdpar_input_error",
    data = list(received = anchor, p = p)
  )
}

#' Internal: resolve the anchor argument to a numeric value
#'
#' Parses the user-supplied \code{anchor} argument of \code{gdpar} and
#' returns the numeric anchor value on the linear-predictor scale.
#'
#' @param anchor Either a numeric scalar or one of
#'   \code{"prior_mean"}, \code{"empirical_y"}.
#' @param family A \code{gdpar_family} object.
#' @param y Numeric vector of outcomes.
#' @param prior A \code{gdpar_prior} object.
#' @param verbose Logical; whether to inform about the chosen anchor.
#' @return Numeric scalar.
#' @keywords internal
#' @noRd
resolve_anchor <- function(anchor, family, y, prior, verbose) {
  if (is.numeric(anchor) && length(anchor) == 1L && is.finite(anchor)) {
    return(as.double(anchor))
  }
  if (is.character(anchor) && length(anchor) == 1L) {
    if (anchor == "prior_mean") {
      return(0.0)
    }
    if (anchor == "empirical_y") {
      yb <- mean(as.double(y))
      eta_anchor <- tryCatch(
        family$linkfun(yb),
        error = function(e) {
          gdpar_abort(
            sprintf(
              "Failed to apply the link function to the outcome mean: %s",
              conditionMessage(e)
            ),
            class = "gdpar_input_error"
          )
        }
      )
      if (verbose) {
        gdpar_inform(
          sprintf(
            "Using empirical anchor: link(%s) of the outcome mean (%g).",
            family$link, eta_anchor
          ),
          class = "gdpar_anchor_message"
        )
      }
      return(eta_anchor)
    }
  }
  gdpar_abort(
    paste(
      "Argument 'anchor' must be a numeric scalar, 'prior_mean'",
      "or 'empirical_y'."
    ),
    class = "gdpar_input_error",
    data = list(received = anchor)
  )
}

#' Internal: collect convergence diagnostics from a cmdstanr fit
#'
#' Runs \code{posterior::summarise_draws} on the relevant parameters,
#' extracts R-hat and effective sample size summaries, and counts
#' divergent transitions and tree-depth saturations from the sampler
#' diagnostics. Issues warnings when thresholds are violated.
#'
#' @param fit A \code{cmdstanr} fit object.
#' @param verbose Logical; whether to issue warnings.
#' @return An object of class \code{gdpar_diagnostics}.
#' @keywords internal
#' @noRd
compute_diagnostics <- function(fit, verbose = TRUE) {
  require_suggested("posterior", "summarize posterior draws")
  draws <- fit$draws()
  vars <- posterior::variables(draws)
  ignore_pat <- "^(eta|log_lik|y_pred|theta_i|a_coef|b_coef|a_raw|b_raw|W_raw)"
  keep_vars <- vars[!grepl(ignore_pat, vars)]
  if (length(keep_vars) == 0L) {
    keep_vars <- intersect(vars, "theta_ref")
    if (length(keep_vars) == 0L) keep_vars <- vars[1L]
  }
  summary <- do.call(
    posterior::summarise_draws,
    c(list(posterior::subset_draws(draws, variable = keep_vars)),
      posterior::default_convergence_measures())
  )
  rhat_vec <- summary$rhat
  ess_bulk_vec <- summary$ess_bulk
  ess_tail_vec <- summary$ess_tail
  rhat_max <- if (length(rhat_vec) > 0L) {
    max(rhat_vec, na.rm = TRUE)
  } else {
    NA_real_
  }
  ess_bulk_min <- if (length(ess_bulk_vec) > 0L) {
    min(ess_bulk_vec, na.rm = TRUE)
  } else {
    NA_real_
  }
  ess_tail_min <- if (length(ess_tail_vec) > 0L) {
    min(ess_tail_vec, na.rm = TRUE)
  } else {
    NA_real_
  }
  divergent_count <- 0L
  treedepth_saturated <- 0L
  efmi_min <- NA_real_
  diag_summary <- tryCatch(fit$diagnostic_summary(),
                            error = function(e) NULL)
  if (!is.null(diag_summary)) {
    if (!is.null(diag_summary$num_divergent)) {
      divergent_count <- as.integer(sum(diag_summary$num_divergent))
    }
    if (!is.null(diag_summary$num_max_treedepth)) {
      treedepth_saturated <- as.integer(sum(diag_summary$num_max_treedepth))
    }
    if (!is.null(diag_summary$ebfmi)) {
      efmi_min <- min(diag_summary$ebfmi, na.rm = TRUE)
    }
  }
  total_iter <- tryCatch(
    as.integer(posterior::niterations(draws)) *
      as.integer(posterior::nchains(draws)),
    error = function(e) NA_integer_
  )
  if (length(total_iter) != 1L || !isTRUE(is.finite(total_iter))) {
    total_iter <- NA_integer_
  }
  rhat_threshold <- 1.01
  ess_threshold <- 400
  div_rel <- if (isTRUE(is.finite(total_iter)) && isTRUE(total_iter > 0L)) {
    divergent_count / total_iter
  } else {
    NA_real_
  }
  tdepth_rel <- if (isTRUE(is.finite(total_iter)) && isTRUE(total_iter > 0L)) {
    treedepth_saturated / total_iter
  } else {
    NA_real_
  }
  converged <- isTRUE(rhat_max <= rhat_threshold) &&
    isTRUE(ess_bulk_min >= ess_threshold) &&
    isTRUE(ess_tail_min >= ess_threshold) &&
    isTRUE(div_rel <= 0.001) &&
    isTRUE(tdepth_rel <= 0.01)
  messages <- character(0)
  if (verbose && !converged) {
    if (isTRUE(rhat_max > rhat_threshold)) {
      messages <- c(messages,
                    sprintf("Maximum R-hat = %.3f exceeds %.2f.",
                            rhat_max, rhat_threshold))
    }
    if (isTRUE(ess_bulk_min < ess_threshold)) {
      messages <- c(messages,
                    sprintf("Minimum bulk ESS = %.0f below %.0f.",
                            ess_bulk_min, ess_threshold))
    }
    if (isTRUE(ess_tail_min < ess_threshold)) {
      messages <- c(messages,
                    sprintf("Minimum tail ESS = %.0f below %.0f.",
                            ess_tail_min, ess_threshold))
    }
    if (isTRUE(div_rel > 0.001)) {
      messages <- c(messages,
                    sprintf("Divergent transitions: %d (%.2f%% of total).",
                            divergent_count, 100 * div_rel))
    }
    if (isTRUE(tdepth_rel > 0.01)) {
      messages <- c(messages,
                    sprintf("Tree-depth saturations: %d (%.2f%% of total).",
                            treedepth_saturated, 100 * tdepth_rel))
    }
    for (m in messages) {
      gdpar_warn(m, class = "gdpar_diagnostic_warning")
    }
  }
  obj <- list(
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergent_count = divergent_count,
    divergent_relative = div_rel,
    treedepth_saturated = treedepth_saturated,
    treedepth_relative = tdepth_rel,
    efmi_min = efmi_min,
    converged = converged,
    summary = summary,
    messages = messages
  )
  class(obj) <- c("gdpar_diagnostics", "list")
  obj
}

#' Extract diagnostics from a fitted gdpar model
#'
#' Convenience accessor for the convergence diagnostics computed at
#' fit time and stored in a \code{gdpar_fit} object.
#'
#' @param fit An object of class \code{gdpar_fit}.
#'
#' @return An object of class \code{gdpar_diagnostics} with R-hat,
#'   effective sample size, divergent-transition and tree-depth
#'   summaries, and a logical \code{converged} flag.
#'
#' @details
#' The diagnostics are computed once at fit time and stored. Calling
#' this function does not re-run any computation.
#'
#' @section Dependencies:
#' Diagnostics are computed with \pkg{posterior} at fit time.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   df <- data.frame(x1 = rnorm(100), y = rnorm(100))
#'   fit <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df,
#'                iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   diagnostics(fit)
#' }
#' }
#'
#' @seealso \code{\link{gdpar}}
#' @export
diagnostics <- function(fit) {
  assert_inherits(fit, "gdpar_fit", "fit")
  fit$diagnostics
}

#' Internal: resolve the user-facing 'group' argument to grouping info
#'
#' Accepts a one-sided formula with a single variable name (for example
#' \code{~ species}) and resolves it against the user-supplied
#' \code{data} into the integer vector consumed by the Stan template.
#' Returns NULL when the user did not supply a grouping. Emits an
#' informational warning when at least one group has fewer than five
#' observations, since the hierarchical anchor is then dominated by
#' shrinkage and the per-group posterior may not be informative.
#'
#' @param group NULL or a one-sided formula identifying a single
#'   variable in \code{data}.
#' @param data User-supplied data frame.
#' @param n Sample size, used for length validation.
#' @param verbose Logical scalar; when TRUE, emit the small-group
#'   informational warning.
#' @return NULL or a list with components \code{group_id} (integer
#'   vector of length n with values in 1..J_groups), \code{var_name}
#'   (character scalar with the user-facing variable name) and
#'   \code{levels} (character vector of the original group levels).
#' @keywords internal
#' @noRd
.resolve_group_argument <- function(group, data, n, verbose = TRUE) {
  if (is.null(group)) return(NULL)
  if (!inherits(group, "formula") || length(group) != 2L) {
    gdpar_abort(
      "Argument 'group' must be a one-sided formula such as ~ species.",
      class = "gdpar_input_error"
    )
  }
  vars <- all.vars(group)
  if (length(vars) != 1L) {
    gdpar_abort(
      paste(
        "Argument 'group' must reference exactly one variable from",
        "'data'; received:",
        paste(vars, collapse = ", "), "."
      ),
      class = "gdpar_input_error",
      data = list(variables_found = vars)
    )
  }
  var_name <- vars[[1L]]
  if (!var_name %in% colnames(data)) {
    gdpar_abort(
      sprintf(
        "Grouping variable '%s' not found in 'data'.", var_name
      ),
      class = "gdpar_input_error",
      data = list(missing_variable = var_name)
    )
  }
  raw <- data[[var_name]]
  if (length(raw) != n) {
    gdpar_abort(
      sprintf(
        paste(
          "Grouping variable '%s' has length %d but the outcome has",
          "length %d; they must coincide."
        ),
        var_name, length(raw), n
      ),
      class = "gdpar_input_error"
    )
  }
  if (any(is.na(raw))) {
    gdpar_abort(
      sprintf(
        "Grouping variable '%s' contains %d NA value(s); please remove or impute them before fitting.",
        var_name, sum(is.na(raw))
      ),
      class = "gdpar_input_error"
    )
  }
  fac <- as.factor(raw)
  J_groups <- nlevels(fac)
  if (J_groups < 1L) {
    gdpar_abort(
      sprintf(
        "Grouping variable '%s' has zero levels.", var_name
      ),
      class = "gdpar_input_error"
    )
  }
  group_id <- as.integer(fac)
  counts <- as.integer(table(group_id))
  if (isTRUE(verbose) && length(counts) > 0L && min(counts) < 5L) {
    gdpar_warn(
      sprintf(
        paste(
          "Grouping variable '%s' has at least one level with fewer",
          "than 5 observations (min n per group = %d, levels = %d).",
          "The per-group anchor theta_ref[g] will be dominated by",
          "shrinkage toward mu_theta_ref; interpret per-group",
          "summaries with caution."
        ),
        var_name, min(counts), J_groups
      ),
      class = "gdpar_grouping_warning"
    )
  }
  list(
    group_id = group_id,
    var_name = var_name,
    levels   = levels(fac)
  )
}

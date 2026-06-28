#' Fit an AMM canonical model via Empirical Bayes (EB)
#'
#' Path 1 Empirical-Bayes counterpart of \code{\link{gdpar}}: estimates
#' the population reference \eqn{\theta_{ref}} by maximizing the
#' marginal likelihood (Type II ML) and samples the lower-level
#' parameters \eqn{\xi = (a, b, W, \sigma_*, \phi)} from the
#' conditional posterior given \eqn{\widehat\theta_{ref}^{EB}}. See
#' \code{vignette("v07_eb_vs_fb", package = "gdpar")} for the
#' Fully-Bayes versus Empirical-Bayes asymptotic comparison
#' (\eqn{p = 1}, \eqn{K = 1}) and
#' \code{vignette("v07b_eb_multivariate", package = "gdpar")} for the
#' multivariate extension that motivates the family of dedicated
#' templates introduced in this sub-phase.
#'
#' @section Scope of Sub-phase 8.6.B:
#' This version implements the base regime \eqn{K = 1}, \eqn{p = 1}
#' (single distributional slot per group, scalar
#' \eqn{\theta_{ref}}). Multivariate \eqn{p > 1} and multi-slot
#' \eqn{K > 1} are explicitly rejected by the input guard with
#' \code{gdpar_unsupported_feature_error}; their habilitation is
#' canonized in Charter Sub-phases 8.6.C and 8.6.D and exercises the
#' \code{K_slots} / \code{p_dim} fields declared in the two dedicated
#' templates \code{inst/stan/amm_eb_marginal.stan} and
#' \code{inst/stan/amm_eb_conditional.stan}.
#'
#' @param formula A two-sided formula \code{y ~ ...}; same semantics as
#'   \code{\link{gdpar}}'s \code{formula} argument.
#' @param family A \code{gdpar_family} object. Sub-phase 8.6.B
#'   exclusively supports the K = 1 univariate families whose canonical
#'   Stan dispatch lives in \code{inst/stan/amm_main.stan}, namely the
#'   four families with \code{stan_id} in \code{c(1, 2, 3, 4)}
#'   (Gaussian, Poisson, neg-binomial-2, Bernoulli). Canonical
#'   K >= 2 families (Beta, Gamma, Student-t, Tweedie,
#'   lognormal-loc-scale) and mixtures (ZIP, ZINB, Hurdle-Poisson,
#'   Hurdle-NB) are rejected with
#'   \code{gdpar_unsupported_feature_error}; their EB habilitation is
#'   wired into Sub-phase 8.6.C jointly with the K > 1 relax.
#' @param amm An \code{amm_spec} object with \code{amm$p == 1L}.
#'   Multivariate specifications are rejected.
#' @param W Optional \code{W_basis} object. Polynomial and B-spline
#'   modulating bases are supported identically to \code{\link{gdpar}}.
#' @param data Data frame containing the variables referenced by
#'   \code{formula} and \code{amm}.
#' @param prior Optional \code{gdpar_prior} object. When \code{NULL},
#'   the package defaults are used.
#' @param anchor Either a numeric scalar or one of
#'   \code{"prior_mean"} (default) and \code{"empirical_y"}.
#' @param skip_id_check Logical scalar; identical semantics to
#'   \code{\link{gdpar}}.
#' @param chains Integer scalar; number of HMC chains for Step (iii)
#'   conditional sampling. Defaults to 4.
#' @param iter_warmup,iter_sampling Integer scalars; HMC warmup and
#'   sampling iterations per chain. Defaults to 1000.
#' @param adapt_delta,max_treedepth,refresh,verbose,seed Identical
#'   semantics to \code{\link{gdpar}}.
#' @param group Optional one-sided formula identifying the grouping
#'   variable in \code{data}. Same semantics as \code{\link{gdpar}}.
#' @param parametrization Character scalar selecting the CP/NCP
#'   sampling parametrization for the additive and modulating
#'   components in Step (iii). One of \code{"auto"} (default; runs the
#'   pre-flight diagnostic), \code{"ncp"}, \code{"cp"}.
#' @param id_check_rigor Character scalar, one of \code{"full"} or
#'   \code{"fast"}.
#' @param eb_correction Logical scalar; when \code{TRUE} (default),
#'   apply the scalar Proposition 7B coverage-discrepancy inflation
#'   factor to the conditional credible intervals (see v07 Section 6 and
#'   v07b Section 5.2). When \code{FALSE}, the credible intervals are
#'   nominal and a \code{gdpar_diagnostic_warning} is issued
#'   advising of the expected \eqn{O(n^{-1})} under-cover.
#' @param laplace_control Named list controlling the Step (i) Laplace
#'   approximation and the anti-fragility strategy of Charter Section 2.8.
#'   Recognized entries (all optional, with documented defaults):
#'   \itemize{
#'     \item \code{multi_start_M} (integer; default 5): number of
#'       independent random inits.
#'     \item \code{kappa_threshold} (numeric; default 1e10): maximum
#'       condition number of the marginal Hessian (after adaptive ridge)
#'       before \code{gdpar_eb_numerical_error} is raised.
#'     \item \code{ridge_init} (numeric; default 1e-6): initial L-M
#'       ridge value used by the adaptive perturbation helper.
#'     \item \code{epsilon_lm} (numeric; default
#'       \code{sqrt(.Machine$double.eps)}, approximately 1.5e-8): the
#'       adaptive Levenberg-Marquardt ridge triggers when either the
#'       posterior covariance has a non-positive or non-finite
#'       eigenvalue, or its determinant is strictly smaller than this
#'       threshold (canon \emph{|det(H)| < epsilon_LM} of Charter
#'       Section 2.8). Canonized in Sub-bloque 9.3.b (Sesion B9.2,
#'       2026-05-27).
#'     \item \code{ridge_max_iter} (integer; default 10): maximum
#'       iterations of the adaptive geometric ridge loop. Each
#'       iteration multiplies the current lambda by
#'       \code{ridge_grow_factor}; the loop terminates with status
#'       \code{"converged"} once the post-ridge condition number is at
#'       or below \code{kappa_threshold}, or with status
#'       \code{"exhausted"} after \code{ridge_max_iter} attempts.
#'     \item \code{ridge_grow_factor} (numeric; default 10.0,
#'       must be > 1): geometric growth factor for the L-M ridge across
#'       iterations.
#'     \item \code{laplace_draws} (integer; default 1000): number of
#'       Gaussian draws produced per Laplace call (used to estimate
#'       the marginal covariance of \eqn{\theta_{ref}}).
#'     \item \code{optim_algorithm} (character; default
#'       \code{"lbfgs"}): algorithm forwarded to
#'       \code{cmdstanr::optimize()}.
#'   }
#' @param ... Additional arguments forwarded to the underlying HMC
#'   sampler of Step (iii).
#'
#' @return An object of class \code{gdpar_eb_fit} with components:
#'   \itemize{
#'     \item \code{theta_ref_hat}: numeric vector of length
#'       \code{J_groups} with the EB point estimates.
#'     \item \code{theta_ref_se}: numeric vector of length
#'       \code{J_groups} with the marginal standard errors derived
#'       from the Laplace covariance.
#'     \item \code{conditional_fit}: the underlying \code{cmdstanr}
#'       fit object of Step (iii).
#'     \item \code{correction_applied}: logical scalar.
#'     \item \code{eb_correction_constant}: numeric scalar; the
#'       scalar Proposition 7B inflation constant when
#'       \code{eb_correction = TRUE}, \code{NA_real_} otherwise.
#'     \item \code{diagnostics_numerical}: list with \code{kappa},
#'       \code{lm_perturbation}, \code{lm_n_iter}, \code{lm_status}
#'       (one of \code{"not_needed"}, \code{"converged"},
#'       \code{"exhausted"}), \code{kappa_post_ridge},
#'       \code{multi_start_dispersion},
#'       \code{marginal_log_lik_history}. For Path C (K > 1 and p > 1)
#'       the slot-vectorized counterparts replace the scalars:
#'       \code{kappa_per_slot}, \code{lm_lambda_per_slot},
#'       \code{lm_n_iter_per_slot}, \code{lm_status_per_slot}.
#'     \item \code{diagnostics}: \code{gdpar_diagnostics} object from
#'       the conditional HMC fit (same shape as in \code{gdpar_fit}).
#'     \item \code{amm}, \code{family}, \code{prior}, \code{design},
#'       \code{anchor}, \code{stan_data}, \code{group_info},
#'       \code{identifiability_report}, \code{parametrization},
#'       \code{call}, \code{path} (\code{"eb"}).
#'   }
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} (Suggests) for Laplace approximation
#' (\code{cmdstanr::laplace()}) in Step (i) and HMC sampling in Step
#' (iii); \pkg{posterior} (Imports) for diagnostics. The Laplace
#' approximation requires the \code{laplace()} method of
#' \code{cmdstanr}, available since cmdstanr 0.7.0.
#'
#' @section References:
#' Carlin, B. P., and Gelfand, A. E. (1990). Approaches for empirical
#' Bayes confidence intervals. JASA 85(409), 105--114.
#'
#' Petrone, S., Rousseau, J., and Scricciolo, C. (2014). Bayes and
#' empirical Bayes: do they merge? Biometrika 101(2), 285--302.
#'
#' Rousseau, J., and Szabo, B. (2017). Asymptotic behaviour of the
#' empirical Bayes posteriors associated to maximum marginal
#' likelihood estimator. Annals of Statistics 45(2), 833--865.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   set.seed(NULL)
#'   n <- 200
#'   df <- data.frame(
#'     x1 = stats::rnorm(n),
#'     x2 = stats::rnorm(n)
#'   )
#'   df$y <- with(df, 1 + 0.5 * x1 - 0.3 * x2 + stats::rnorm(n, sd = 0.5))
#'   spec <- amm_spec(a = ~ x1 + x2)
#'   fit_eb <- gdpar_eb(
#'     formula       = y ~ x1 + x2,
#'     family        = gdpar_family("gaussian"),
#'     amm           = spec,
#'     data          = df,
#'     iter_warmup   = 200,
#'     iter_sampling = 200,
#'     chains        = 2
#'   )
#'   print(fit_eb)
#' }
#' }
#'
#' @seealso \code{\link{gdpar}}, \code{\link{amm_spec}},
#'   \code{\link{gdpar_family}}, \code{\link{gdpar_prior}}.
#' @export
gdpar_eb <- function(formula,
                     family = gdpar_family("gaussian"),
                     amm = amm_spec(),
                     W = NULL,
                     data,
                     prior = NULL,
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
                     id_check_rigor = c("full", "fast"),
                     eb_correction = TRUE,
                     laplace_control = list(),
                     ...) {
  parametrization <- match.arg(parametrization)
  id_check_rigor  <- match.arg(id_check_rigor)
  call <- match.call()

  .gdpar_eb_validate_inputs(
    formula = formula, family = family, amm = amm, data = data,
    eb_correction = eb_correction, laplace_control = laplace_control
  )

  if (is.null(prior)) prior <- gdpar_prior()
  assert_inherits(prior, "gdpar_prior", "prior")

  lc <- .gdpar_eb_resolve_laplace_control(laplace_control)

  # Sub-phase 8.6.C (decision D34) input dispatch: detect the three
  # K-input patterns mirroring gdpar() (formula = gdpar_formula_set;
  # amm = named list of amm_spec; legacy formula with a()/b()/W() in
  # its RHS). When any of them fires we build amm_list_canonical, then
  # either route to the Path B orchestrator .gdpar_eb_run_K() (when
  # the resolved K > 1) or unwrap to a single amm_spec and fall through
  # to the K = 1 pipeline below (Path A p > 1 still flows through that
  # path because assemble_stan_data() dispatches to its multivariate
  # branch internally when amm$p > 1).
  .formula_set_input <- inherits(formula, "gdpar_formula_set")
  .amm_list_input <- is.list(amm) && !inherits(amm, "amm_spec") &&
                      !is.null(names(amm))
  .classic_with_amm_calls <- !.formula_set_input &&
                              inherits(formula, "formula") &&
                              length(formula) == 3L &&
                              .gdpar_rhs_has_amm_calls(formula)
  .family_is_named_list <- is.list(family) &&
                            !inherits(family, "gdpar_family") &&
                            !inherits(family, "gdpar_family_multi") &&
                            !is.null(names(family))
  if (.formula_set_input || .amm_list_input || .classic_with_amm_calls) {
    K_inputs <- .gdpar_eb_resolve_K_inputs(
      formula = formula, amm = amm, W = W, family = family,
      formula_set_input = .formula_set_input,
      amm_list_input = .amm_list_input,
      classic_with_amm_calls = .classic_with_amm_calls,
      family_is_named_list = .family_is_named_list
    )
    if (K_inputs$K > 1L) {
      # Sub-phase 8.6.D (decision D38'' = (h) canonized Session 13b
      # 2026-05-25): detect Path C by inspecting whether any slot in
      # amm_list_canonical carries p > 1; if so, route to
      # .gdpar_eb_run_KxP() instead of .gdpar_eb_run_K().
      .any_slot_p_gt1 <- any(vapply(
        K_inputs$amm_list_canonical,
        function(a) isTRUE(!is.null(a$p) && a$p > 1L),
        logical(1L)
      ))
      if (.any_slot_p_gt1) {
        return(.gdpar_eb_run_KxP(
          amm_list_canonical = K_inputs$amm_list_canonical,
          family             = K_inputs$family_promoted,
          data               = data,
          prior              = prior,
          anchor             = anchor,
          outcome_name       = K_inputs$outcome_name,
          formula_env        = K_inputs$formula_env,
          family_id_k_vector = K_inputs$family_id_k_vector,
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
          id_check_rigor     = id_check_rigor,
          eb_correction      = eb_correction,
          laplace_control    = lc,
          call               = call,
          ...
        ))
      }
      return(.gdpar_eb_run_K(
        amm_list_canonical = K_inputs$amm_list_canonical,
        family             = K_inputs$family_promoted,
        data               = data,
        prior              = prior,
        anchor             = anchor,
        outcome_name       = K_inputs$outcome_name,
        formula_env        = K_inputs$formula_env,
        family_id_k_vector = K_inputs$family_id_k_vector,
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
        id_check_rigor     = id_check_rigor,
        eb_correction      = eb_correction,
        laplace_control    = lc,
        call               = call,
        ...
      ))
    }
    amm <- K_inputs$amm_list_canonical[[1L]]
    family <- K_inputs$family_promoted
    .vars_for_rhs <- union(
      if (!is.null(amm$a)) all.vars(amm$a) else character(0L),
      if (!is.null(amm$b)) all.vars(amm$b) else character(0L)
    )
    .rhs_str <- if (length(.vars_for_rhs) == 0L) "1" else
                paste(.vars_for_rhs, collapse = " + ")
    formula <- stats::as.formula(
      paste(K_inputs$outcome_name, "~", .rhs_str),
      env = K_inputs$formula_env
    )
  }

  # Resolved (K, p) for the K = 1 pipeline (Path 8.6.B when p = 1,
  # Path A when amm$p > 1). The dispatcher of Stan templates and the
  # K_slots / p_dim fields injected into stan_data both read from
  # these resolved scalars.
  p_resolved <- if (inherits(amm, "amm_spec") && !is.null(amm$p)) {
    as.integer(amm$p)
  } else {
    1L
  }
  K_resolved <- 1L
  .gdpar_eb_check_stan_id_for_path(family, K_resolved, p_resolved)

  outcome_name <- as.character(formula[[2]])
  if (!outcome_name %in% colnames(data)) {
    gdpar_abort(
      sprintf("Outcome variable '%s' not found in data.", outcome_name),
      class = "gdpar_input_error"
    )
  }
  y <- data[[outcome_name]]
  y_bad <- if (is.numeric(y)) !is.finite(y) else is.na(y)
  if (any(y_bad)) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' contains %d non-finite value(s) (NA, NaN or Inf).",
        outcome_name, sum(y_bad)
      ),
      class = "gdpar_input_error"
    )
  }
  rhs <- formula[c(1L, 3L)]
  rhs <- stats::update(rhs, ~ . + 0)
  if (!is.null(amm$W)) {
    amm$W <- materialize_W_basis(amm$W, p = p_resolved)
  }
  design <- build_amm_design(amm, data, formula_rhs = rhs)
  anchor_value <- resolve_anchor(anchor, family, y, prior, verbose)

  id_report <- if (skip_id_check) NULL else {
    diag_theta <- if (!is.null(amm$b) && abs(anchor_value) < 1e-8) 1
                  else anchor_value
    rep <- gdpar_check_identifiability(
      amm = amm, data = data,
      theta_ref_init = diag_theta,
      formula_rhs = rhs
    )
    if (!isTRUE(rep$passed)) {
      gdpar_abort(
        paste(
          "Basis-restricted identifiability check failed at the",
          "diagnostic test point. Inspect the returned report or pass",
          "skip_id_check = TRUE to override at your own risk."
        ),
        class = "gdpar_identifiability_error",
        data = list(report = rep)
      )
    }
    rep
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
  stan_data$K_slots <- as.integer(K_resolved)
  stan_data$p_dim   <- as.integer(p_resolved)

  parametrization_resolved <- resolve_parametrization(
    parametrization      = parametrization,
    parametrization_a    = NULL,
    parametrization_W    = NULL,
    prior                = prior,
    stan_data            = stan_data,
    amm                  = amm,
    preflight_seed       = seed,
    verbose              = verbose
  )

  marginal_src <- .gdpar_eb_generate_stan_marginal(
    prior,
    cp_a   = parametrization_resolved$cp_a,
    cp_W   = parametrization_resolved$cp_W,
    K      = K_resolved,
    p      = p_resolved,
    family = family
  )
  marginal_path <- write_stan_to_tempfile(marginal_src)
  require_suggested("cmdstanr",
                    "fit Path 1 Empirical-Bayes models via cmdstanr::laplace()")
  marginal_model <- cmdstanr::cmdstan_model(marginal_path)

  laplace_result <- .gdpar_eb_maximize_marginal(
    model     = marginal_model,
    stan_data = stan_data,
    control   = lc,
    seed      = seed,
    verbose   = verbose
  )

  conditional_src <- .gdpar_eb_generate_stan_conditional(
    prior,
    cp_a   = parametrization_resolved$cp_a,
    cp_W   = parametrization_resolved$cp_W,
    K      = K_resolved,
    p      = p_resolved,
    family = family
  )
  conditional_path  <- write_stan_to_tempfile(conditional_src)
  conditional_model <- cmdstanr::cmdstan_model(conditional_path)

  stan_data_cond <- stan_data
  # Path-aware shape for theta_ref_data (Sub-phase 8.6.C D35 emergent
  # 2026-05-25 Session 12): the conditional template declares
  # array[J_groups] vector[p] theta_ref_data. For K = 1 + p = 1 the
  # flat scalar consumed by 8.6.B preserves bit-exactness; for
  # Path A (K = 1, p > 1) the Laplace helper returns a length
  # J_groups * p numeric vector that must be reshaped to a
  # J_groups by p matrix (cmdstanr's automatic packing accepts both
  # nested lists and matrices for array[J] vector[p]).
  J_groups_loc <- as.integer(stan_data$J_groups)
  theta_hat_loc <- as.numeric(laplace_result$theta_ref_hat)
  if (p_resolved > 1L &&
      length(theta_hat_loc) == J_groups_loc * p_resolved) {
    stan_data_cond$theta_ref_data <- matrix(theta_hat_loc,
                                            nrow = J_groups_loc,
                                            ncol = p_resolved,
                                            byrow = FALSE)
  } else {
    stan_data_cond$theta_ref_data <- theta_hat_loc
  }

  sample_args <- list(
    data            = stan_data_cond,
    chains          = as.integer(chains),
    iter_warmup     = as.integer(iter_warmup),
    iter_sampling   = as.integer(iter_sampling),
    adapt_delta     = adapt_delta,
    max_treedepth   = as.integer(max_treedepth),
    refresh         = as.integer(refresh),
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) sample_args$seed <- as.integer(seed)
  extra <- list(...)
  for (nm in names(extra)) sample_args[[nm]] <- extra[[nm]]
  fit_cond <- do.call(conditional_model$sample, sample_args)
  diagnostics <- compute_diagnostics(fit_cond, verbose = verbose)

  correction <- .gdpar_eb_apply_correction(
    eb_correction  = eb_correction,
    laplace_result = laplace_result,
    stan_data      = stan_data,
    p              = p_resolved,
    verbose        = verbose
  )

  obj <- list(
    theta_ref_hat            = laplace_result$theta_ref_hat,
    theta_ref_se             = laplace_result$theta_ref_se,
    conditional_fit          = fit_cond,
    amm                      = amm,
    family                   = family,
    prior                    = prior,
    design                   = design,
    anchor                   = anchor_value,
    stan_data                = stan_data,
    identifiability_report   = id_report,
    diagnostics              = diagnostics,
    diagnostics_numerical    = laplace_result$diagnostics,
    parametrization          = list(
      cp_a = parametrization_resolved$cp_a,
      cp_W = parametrization_resolved$cp_W,
      meta = parametrization_resolved$meta
    ),
    group_info               = group_info,
    correction_applied       = correction$applied,
    eb_correction_constant   = correction$constant,
    call                     = call,
    path                     = "eb"
  )
  class(obj) <- c("gdpar_eb_fit", "list")
  obj
}

#' Internal: validate the inputs to gdpar_eb() before the heavy machinery
#'
#' Sub-phase 8.6.C (decision D34) + Sub-phase 8.6.D (decisions D36 =
#' alpha + D38 = (b'), Session 13b 2026-05-25): enforces the type
#' discipline of the high-level arguments (\code{formula} is a formula
#' or \code{gdpar_formula_set}; \code{amm} is an \code{amm_spec} or
#' named list of \code{amm_spec}; \code{data} is a data frame; etc.)
#' Sub-phase 8.6.D releases the last guarded combined regime
#' \code{K > 1 AND p > 1}: the resolved (K, p) is routed to the new
#' Path C orchestrator \code{.gdpar_eb_run_KxP()}. The fine-grained
#' per-path checks (supported \code{stan_id} for the resolved (K, p),
#' \code{family$min_K} consistency, etc.) live in the
#' \code{gdpar_eb()} dispatcher after the K-input pattern is resolved
#' so the error messages can speak in terms of the actually resolved
#' Path A / Path B / Path C regime.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_validate_inputs <- function(formula, family, amm, data,
                                       eb_correction, laplace_control) {
  .formula_set_input <- inherits(formula, "gdpar_formula_set")
  if (!.formula_set_input &&
      (!inherits(formula, "formula") || length(formula) != 3L)) {
    gdpar_abort(
      paste(
        "Argument 'formula' must be a two-sided formula 'y ~ ...' or a",
        "'gdpar_formula_set' object."
      ),
      class = "gdpar_input_error"
    )
  }
  .family_is_named_list <- is.list(family) &&
                            !inherits(family, "gdpar_family") &&
                            !inherits(family, "gdpar_family_multi") &&
                            !is.null(names(family)) &&
                            all(nzchar(names(family))) &&
                            anyDuplicated(names(family)) == 0L &&
                            all(vapply(family, inherits, logical(1L),
                                       "gdpar_family"))
  if (!.family_is_named_list &&
      !inherits(family, c("gdpar_family", "gdpar_family_multi"))) {
    gdpar_abort(
      paste(
        "Argument 'family' must be a 'gdpar_family' object, a",
        "'gdpar_family_multi' object (Path A, p > 1), or a named list",
        "of 'gdpar_family' objects (Path B heterogeneous K, sub-phase",
        "8.3.7 pattern)."
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(family))
    )
  }
  .amm_list_input <- is.list(amm) && !inherits(amm, "amm_spec") &&
                      !is.null(names(amm))
  if (!.amm_list_input && !inherits(amm, "amm_spec")) {
    gdpar_abort(
      paste(
        "Argument 'amm' must be an 'amm_spec' object or a named list of",
        "'amm_spec' objects (Path B K > 1)."
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(amm))
    )
  }
  # Sub-phase 8.6.D (Session 13b 2026-05-25): the last K > 1 + p > 1
  # guard is released; Path C is routed to .gdpar_eb_run_KxP() in the
  # dispatcher below. Per-path supported sets (stan_id, min_K,
  # heterogeneous slot patterns) are checked downstream.
  if (!is.data.frame(data)) {
    gdpar_abort(
      "Argument 'data' must be a data frame.",
      class = "gdpar_input_error"
    )
  }
  if (!is.logical(eb_correction) || length(eb_correction) != 1L ||
      is.na(eb_correction)) {
    gdpar_abort(
      "Argument 'eb_correction' must be a single non-NA logical value.",
      class = "gdpar_input_error"
    )
  }
  if (!is.list(laplace_control)) {
    gdpar_abort(
      "Argument 'laplace_control' must be a (possibly empty) named list.",
      class = "gdpar_input_error"
    )
  }
  invisible(NULL)
}

#' Internal: enforce the per-path supported stan_id set
#'
#' Sub-phase 8.6.C (decision D33 relax) + Sub-phase 8.6.D (decisions
#' D38'' = (h) + D40' canonized Session 13b 2026-05-25): the
#' EB-supported \code{stan_id} set depends on the resolved (K, p)
#' regime:
#' \itemize{
#'   \item (K = 1, p = 1): {1, 2, 3, 4} (Gaussian, Poisson, NB,
#'     Bernoulli) -- amm_eb_marginal.stan body inherits amm_main.stan.
#'   \item (K = 1, p > 1): {1, 2, 3, 4} (homogeneous coord-wise) --
#'     amm_eb_marginal_multi.stan body inherits amm_distrib_multi.stan.
#'   \item (K > 1, p = 1): {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13} --
#'     amm_eb_marginal_K.stan body inherits amm_distrib_K.stan with the
#'     full distributional-regression dispatch.
#'   \item (K > 1, p > 1): {1, 3} (Gaussian K=2, NB K=2) --
#'     amm_eb_marginal_KxP.stan with coord-wise factorization
#'     composed with Path B K-slot multi-parametric likelihood
#'     (canonical composition of Theorem 7C* + Theorem 7A* of v07b).
#'     The remaining Path B set {5, 6, 7, 8, 9, 10, 11, 12, 13}
#'     (Beta, Gamma, Lognormal_loc_scale, Student-t, Tweedie, ZIP,
#'     ZINB, Hurdle-Poisson, Hurdle-NB) is deferred to a later
#'     iteration of 8.6.D with the explicit numerical caveat of
#'     opening Section 6.1 (HMC condicional bajo plug-in EB cerca del
#'     borde de soporte logit/log links + warmup corto).
#' }
#'
#' Under Path C the dispatcher iterates the check across the K slots
#' before assembling the family_id_k_vector data field.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_check_stan_id_for_path <- function(family, K, p) {
  if (is.null(family$stan_id)) return(invisible(NULL))
  stan_id <- as.integer(family$stan_id)
  supported <- if (K == 1L) {
    c(1L, 2L, 3L, 4L)
  } else if (p > 1L) {
    # Path C (Sub-phase 8.6.D first iteration, decision D40'):
    # K-slot Path B families restricted to {1, 3} = Gaussian K=2 + NB
    # K=2. Other Path B families deferred per opening Section 6.1.
    c(1L, 3L)
  } else {
    # Path B: full K-individual family set.
    c(1L, 3L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L)
  }
  if (!isTRUE(stan_id %in% supported)) {
    gdpar_abort(
      sprintf(
        paste(
          "Family '%s' (stan_id = %d) is not wired in the EB template",
          "for the resolved regime (K = %d, p = %d). Supported stan_ids",
          "in this regime: %s."
        ),
        family$name, stan_id, K, p,
        paste(supported, collapse = ", ")
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(family = family$name, stan_id = stan_id,
                  K = K, p = p, supported = supported)
    )
  }
  invisible(NULL)
}

#' Internal: merge user-supplied laplace_control with documented defaults
#'
#' Recognized entries are \code{multi_start_M}, \code{kappa_threshold},
#' \code{ridge_init}, \code{laplace_draws}, \code{optim_algorithm}.
#' Unknown entries trigger a soft warning rather than an error so users
#' can experiment with forward-compatible options without rewriting
#' their call.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_resolve_laplace_control <- function(user) {
  defaults <- list(
    multi_start_M     = 5L,
    kappa_threshold   = 1e10,
    ridge_init        = 1e-6,
    laplace_draws     = 1000L,
    optim_algorithm   = "lbfgs",
    epsilon_lm        = sqrt(.Machine$double.eps),
    ridge_max_iter    = 10L,
    ridge_grow_factor = 10.0
  )
  if (length(user) == 0L) return(defaults)
  if (is.null(names(user)) || any(!nzchar(names(user)))) {
    gdpar_abort(
      "Every entry of 'laplace_control' must be named.",
      class = "gdpar_input_error"
    )
  }
  unknown <- setdiff(names(user), names(defaults))
  if (length(unknown) > 0L) {
    gdpar_warn(
      sprintf(
        paste(
          "Unknown 'laplace_control' entries ignored: %s. Recognized:",
          "multi_start_M, kappa_threshold, ridge_init, laplace_draws,",
          "optim_algorithm, epsilon_lm, ridge_max_iter,",
          "ridge_grow_factor."
        ),
        paste(sQuote(unknown), collapse = ", ")
      ),
      class = "gdpar_diagnostic_warning"
    )
  }
  out <- defaults
  for (nm in intersect(names(user), names(defaults))) {
    out[[nm]] <- user[[nm]]
  }
  out$multi_start_M     <- as.integer(out$multi_start_M)
  out$laplace_draws     <- as.integer(out$laplace_draws)
  out$kappa_threshold   <- as.double(out$kappa_threshold)
  out$ridge_init        <- as.double(out$ridge_init)
  out$epsilon_lm        <- as.double(out$epsilon_lm)
  out$ridge_max_iter    <- as.integer(out$ridge_max_iter)
  out$ridge_grow_factor <- as.double(out$ridge_grow_factor)
  if (!isTRUE(out$multi_start_M >= 1L)) {
    gdpar_abort(
      "'laplace_control$multi_start_M' must be a positive integer.",
      class = "gdpar_input_error"
    )
  }
  if (!isTRUE(out$kappa_threshold > 0)) {
    gdpar_abort(
      "'laplace_control$kappa_threshold' must be strictly positive.",
      class = "gdpar_input_error"
    )
  }
  if (!isTRUE(out$epsilon_lm > 0)) {
    gdpar_abort(
      "'laplace_control$epsilon_lm' must be strictly positive.",
      class = "gdpar_input_error"
    )
  }
  if (!isTRUE(out$ridge_max_iter >= 1L)) {
    gdpar_abort(
      "'laplace_control$ridge_max_iter' must be a positive integer.",
      class = "gdpar_input_error"
    )
  }
  if (!isTRUE(out$ridge_grow_factor > 1)) {
    gdpar_abort(
      "'laplace_control$ridge_grow_factor' must be strictly greater than 1.",
      class = "gdpar_input_error"
    )
  }
  out
}

#' Internal: adaptive Levenberg-Marquardt ridge for the empirical Laplace covariance
#'
#' Canonizes component 2 of the four-component anti-fragility strategy
#' of \code{CHARTER_SUBFASE_8_6.md} Section 2.8, extending the
#' single-step ridge of Sub-phase 8.6.B into an adaptive geometric loop
#' canonized in Sub-bloque 9.3.b of Bloque 9 (Sesion B9.2, 2026-05-27).
#'
#' Triggers on either: (a) any eigenvalue non-finite or non-positive,
#' or (b) |det(cov)| < control$epsilon_lm. When triggered, an effective
#' ridge \eqn{\lambda_{\mathrm{eff}} = \max(\lambda, 10^{-3} |\bar{\mathrm{diag}}|)}
#' is applied; if the resulting matrix is not positive-definite or its
#' condition number remains above \code{control$kappa_threshold},
#' \eqn{\lambda} grows by \code{control$ridge_grow_factor} and the loop
#' iterates up to \code{control$ridge_max_iter} times.
#'
#' @param cov Numeric square symmetric matrix; the empirical posterior
#'   covariance returned by \code{cmdstanr::laplace()} (per-slot block
#'   for Path C, full theta_ref covariance for Path A/B).
#' @param control Resolved \code{laplace_control} list.
#'
#' @return List with fields:
#'   \describe{
#'     \item{\code{cov_perturbed}}{Possibly-ridged covariance (or
#'       original when status is \code{"not_needed"}).}
#'     \item{\code{lambda_used}}{Final effective ridge (0 when
#'       not_needed).}
#'     \item{\code{n_iter}}{Iterations performed (0 when not_needed).}
#'     \item{\code{kappa_post}}{Condition number after ridge (or
#'       original kappa when not_needed; Inf when exhausted).}
#'     \item{\code{status}}{Character scalar in
#'       \code{c("not_needed", "converged", "exhausted")}.}
#'   }
#'
#' @keywords internal
#' @noRd
.gdpar_eb_lm_perturb <- function(cov, control) {
  dim_n <- nrow(cov)
  trace_mean <- max(abs(mean(diag(cov))), 1e-12)
  eigs0 <- tryCatch(
    eigen(cov, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  det_val <- if (all(is.finite(eigs0))) prod(eigs0) else NA_real_
  needs_ridge <- (any(!is.finite(eigs0)) || any(eigs0 <= 0) ||
                   (is.finite(det_val) &&
                    abs(det_val) < control$epsilon_lm))
  if (!needs_ridge) {
    kappa_pre <- max(eigs0) / min(eigs0)
    return(list(
      cov_perturbed = cov,
      lambda_used   = 0,
      n_iter        = 0L,
      kappa_post    = kappa_pre,
      status        = "not_needed"
    ))
  }
  lambda <- control$ridge_init
  cov_try <- cov
  eigs_try <- eigs0
  lambda_eff <- 0
  for (it in seq_len(control$ridge_max_iter)) {
    lambda_eff <- max(lambda, 1e-3 * trace_mean)
    cov_try <- cov + diag(lambda_eff, dim_n, dim_n)
    eigs_try <- tryCatch(
      eigen(cov_try, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    if (all(is.finite(eigs_try)) && min(eigs_try) > 0) {
      kappa_try <- max(eigs_try) / min(eigs_try)
      if (kappa_try <= control$kappa_threshold) {
        return(list(
          cov_perturbed = cov_try,
          lambda_used   = lambda_eff,
          n_iter        = it,
          kappa_post    = kappa_try,
          status        = "converged"
        ))
      }
    }
    lambda <- lambda * control$ridge_grow_factor
  }
  kappa_final <- if (all(is.finite(eigs_try)) && min(eigs_try) > 0) {
    max(eigs_try) / min(eigs_try)
  } else Inf
  list(
    cov_perturbed = cov_try,
    lambda_used   = lambda_eff,
    n_iter        = control$ridge_max_iter,
    kappa_post    = kappa_final,
    status        = "exhausted"
  )
}

#' Internal: render the EB marginal Stan template for the resolved (K, p)
#'
#' Sub-phase 8.6.C (decision D34, canonized Session 10): dispatch on
#' (K, p) selects one of four EB templates parallel to the FB triplet:
#' \itemize{
#'   \item (K = 1, p = 1): \code{amm_eb_marginal.stan} (Sub-phase 8.6.B,
#'     bit-identical body to \code{amm_main.stan}).
#'   \item (K = 1, p > 1): \code{amm_eb_marginal_multi.stan} (D34
#'     Path A, bit-identical body to \code{amm_distrib_multi.stan}).
#'   \item (K > 1, p = 1): \code{amm_eb_marginal_K.stan} (D34 Path B,
#'     bit-identical body to \code{amm_distrib_K.stan}).
#'   \item (K > 1, p > 1): guarded upstream
#'     (\code{gdpar_unsupported_feature_error}; queued for 8.6.D).
#' }
#'
#' For the K = 1 + p = 1 leaf the helper preserves the byte-identical
#' source of 8.6.B by routing through \code{.gdpar_eb_render_template}.
#' For the multi and K leaves the helper delegates to the FB generators
#' \code{generate_stan_code_multi} / \code{generate_stan_code_K} with
#' the \code{template_name} override added in Sub-phase 8.6.C (decision
#' D34), guaranteeing that EB Path A and Path B share the exact
#' placeholder substitution logic as their FB counterparts.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_generate_stan_marginal <- function(prior, cp_a = FALSE,
                                              cp_W = FALSE,
                                              K = 1L, p = 1L,
                                              family = NULL,
                                              cp_a_per_k = NULL,
                                              cp_a_per_K = NULL) {
  K <- as.integer(K)
  p <- as.integer(p)
  if (K == 1L && p == 1L) {
    return(.gdpar_eb_render_template("amm_eb_marginal.stan",
                                     prior, cp_a, cp_W))
  }
  if (K == 1L && p > 1L) {
    return(generate_stan_code_multi(
      prior, cp_a = cp_a, cp_W = cp_W,
      cp_a_per_k    = cp_a_per_k,
      template_name = "amm_eb_marginal_multi.stan"
    ))
  }
  if (K > 1L && p == 1L) {
    return(generate_stan_code_K(
      prior, cp_a = cp_a, cp_W = cp_W,
      cp_a_per_K    = cp_a_per_K,
      family        = family,
      template_name = "amm_eb_marginal_K.stan"
    ))
  }
  # Sub-phase 8.6.D (Session 13b): Path C K > 1 AND p > 1 via dedicated
  # template (decision D36 = alpha). Placeholders restricted to the
  # prior set: theta_ref, sigma_theta_ref, sigma_a, sigma_b, sigma_y,
  # phi. {{A_SCALE}}, {{A_PRIOR}}, {{W_SCALE}}, {{W_PRIOR}} are absent
  # from the Path C template body (NCP per slot per coord hardcoded;
  # W disabled per D39).
  return(.gdpar_eb_render_template("amm_eb_marginal_KxP.stan",
                                   prior, cp_a, cp_W))
}

#' Internal: render the EB conditional Stan template for the resolved (K, p)
#'
#' Sub-phase 8.6.C (decision D34): companion of
#' \code{.gdpar_eb_generate_stan_marginal} for Step (iii) of the EB
#' workflow. Dispatch table identical to the marginal helper; the
#' conditional templates move \code{theta_ref} (or \code{theta_ref_k})
#' from \code{parameters\{\}} to \code{data\{\}} and drop the anchor
#' priors from \code{model\{\}} but otherwise share the placeholder set
#' with the marginal templates.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_generate_stan_conditional <- function(prior, cp_a = FALSE,
                                                 cp_W = FALSE,
                                                 K = 1L, p = 1L,
                                                 family = NULL,
                                                 cp_a_per_k = NULL,
                                                 cp_a_per_K = NULL) {
  K <- as.integer(K)
  p <- as.integer(p)
  if (K == 1L && p == 1L) {
    return(.gdpar_eb_render_template("amm_eb_conditional.stan",
                                     prior, cp_a, cp_W))
  }
  if (K == 1L && p > 1L) {
    return(generate_stan_code_multi(
      prior, cp_a = cp_a, cp_W = cp_W,
      cp_a_per_k    = cp_a_per_k,
      template_name = "amm_eb_conditional_multi.stan"
    ))
  }
  if (K > 1L && p == 1L) {
    return(generate_stan_code_K(
      prior, cp_a = cp_a, cp_W = cp_W,
      cp_a_per_K    = cp_a_per_K,
      family        = family,
      template_name = "amm_eb_conditional_K.stan"
    ))
  }
  # Sub-phase 8.6.D (Session 13b): Path C K > 1 AND p > 1.
  return(.gdpar_eb_render_template("amm_eb_conditional_KxP.stan",
                                   prior, cp_a, cp_W))
}

#' Internal: shared template renderer for the two EB Stan templates
#'
#' Reproduces the placeholder substitution logic of
#' \code{generate_stan_code()} verbatim, restricted to the EB template
#' family. Aborts with a structured error if a placeholder remains
#' un-substituted, matching the safety net of the FB codegen.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_render_template <- function(template_name, prior, cp_a, cp_W) {
  # Sub-bloque 9.3.a colateral L.iv.A.2 (Bloque 9, Sesion B9.6,
  # 2026-05-27): translate the EB-single legacy template names to their
  # canonical-piece equivalents in inst/stan/_canonical_pieces/. The
  # KxP EB templates (amm_eb_marginal_KxP.stan +
  # amm_eb_conditional_KxP.stan) canonized in Sub-fase 8.6.D do NOT
  # carry inline bspline + apply_W_basis_diff helpers (use_W = 0 per
  # D39); they stay in inst/stan/ root unchanged by the cascade.
  effective_template <- switch(template_name,
    "amm_eb_marginal.stan"     = "amm_canonical_eb_marginal.stan",
    "amm_eb_conditional.stan"  = "amm_canonical_eb_conditional.stan",
    template_name
  )
  if (startsWith(effective_template, "amm_canonical_")) {
    template_path <- system.file("stan", "_canonical_pieces",
                                 effective_template, package = "gdpar")
    if (!nzchar(template_path)) {
      template_path <- file.path("inst", "stan", "_canonical_pieces",
                                 effective_template)
    }
  } else {
    template_path <- system.file("stan", effective_template, package = "gdpar")
    if (!nzchar(template_path)) {
      template_path <- file.path("inst", "stan", effective_template)
    }
  }
  if (!file.exists(template_path)) {
    gdpar_abort(
      sprintf("Stan template file '%s' not found.", effective_template),
      class = "gdpar_internal_error"
    )
  }
  src <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
  # Sub-bloque 9.3.a colateral L.iv.A.2 (Sesion B9.6, 2026-05-27):
  # mirror the helpers substitution that `.gdpar_emit_canonical_stan()`
  # already performs for the FB dispatcher. When the body piece contains
  # the `// {{CANONICAL_HELPERS}}` placeholder, inject the canonical
  # helpers piece pre-substitution. Forward-compatible: pieces without
  # the placeholder (KxP EB) pass through unchanged.
  if (grepl("// {{CANONICAL_HELPERS}}", src, fixed = TRUE)) {
    helpers_path <- system.file("stan", "_canonical_pieces",
                                "amm_canonical_helpers.stan",
                                package = "gdpar")
    if (!nzchar(helpers_path)) {
      helpers_path <- file.path("inst", "stan", "_canonical_pieces",
                                "amm_canonical_helpers.stan")
    }
    if (!file.exists(helpers_path)) {
      gdpar_abort(
        "Canonical helpers piece 'amm_canonical_helpers.stan' not found.",
        class = "gdpar_internal_error"
      )
    }
    helpers_src <- paste(readLines(helpers_path, warn = FALSE),
                         collapse = "\n")
    src <- gsub("// {{CANONICAL_HELPERS}}", helpers_src, src, fixed = TRUE)
  }
  replacements <- list(
    "{{PRIOR_THETA_REF}}"       = prior$theta_ref,
    "{{PRIOR_SIGMA_THETA_REF}}" = prior$sigma_theta_ref,
    "{{PRIOR_SIGMA_A}}"         = prior$sigma_a,
    "{{PRIOR_SIGMA_B}}"         = prior$sigma_b,
    "{{PRIOR_SIGMA_W}}"         = prior$sigma_W,
    "{{PRIOR_SIGMA_Y}}"         = prior$sigma_y,
    "{{PRIOR_PHI}}"             = prior$phi,
    "{{A_SCALE}}"               = if (cp_a) "" else " * sigma_a[1]",
    "{{A_PRIOR}}"               = if (cp_a) "normal(0, sigma_a[1])"
                                         else "normal(0, 1)",
    "{{W_SCALE}}"               = if (cp_W) "" else " * sigma_W[1]",
    "{{W_PRIOR}}"               = if (cp_W) "normal(0, sigma_W[1])"
                                         else "normal(0, 1)"
  )
  for (k in names(replacements)) {
    src <- gsub(k, replacements[[k]], src, fixed = TRUE)
  }
  if (grepl("{{", src, fixed = TRUE)) {
    leftover <- regmatches(
      src, regexpr("\\{\\{[A-Za-z0-9_]+\\}\\}", src)
    )
    gdpar_abort(
      sprintf(
        "Unsubstituted placeholder remains in EB Stan code: %s",
        leftover
      ),
      class = "gdpar_internal_error"
    )
  }
  src
}

#' Internal: Step (i) of the EB workflow with the anti-fragility
#' strategy of Charter Section 2.8
#'
#' Runs \code{cmdstanr::optimize()} + \code{cmdstanr::laplace()} on
#' \code{amm_eb_marginal.stan} with \code{multi_start_M} independent
#' random inits, retains the init with the highest log-marginal
#' approximation, applies the adaptive Levenberg-Marquardt ridge if
#' the Hessian is singular, and reports the diagnostics needed by the
#' \code{gdpar_eb_fit$diagnostics_numerical} slot.
#'
#' @return A list with components \code{theta_ref_hat} (numeric vector
#'   length \code{J_groups}), \code{theta_ref_se} (numeric vector
#'   length \code{J_groups}), \code{theta_ref_cov} (matrix
#'   \code{J_groups x J_groups}), and \code{diagnostics} (named list
#'   with \code{kappa}, \code{lm_perturbation},
#'   \code{multi_start_dispersion}, \code{marginal_log_lik_history}).
#'   Aborts with \code{gdpar_eb_numerical_error} when all inits
#'   produce singular Hessians or with
#'   \code{gdpar_unsupported_feature_error} when the Laplace fallback
#'   recommends FB.
#' @keywords internal
#' @noRd
.gdpar_eb_maximize_marginal <- function(model, stan_data, control,
                                         seed, verbose) {
  M <- control$multi_start_M
  J <- stan_data$J_groups
  history_lp <- rep(NA_real_, M)
  best_idx <- NA_integer_
  best_opt <- NULL
  best_init <- NULL
  # Path-aware init dispatch (Sub-phase 8.6.C, decision D35 emergent
  # 2026-05-25 Session 12): the K = 1 + p = 1 random_init helper is
  # bit-exact preserved; for Path A (p > 1) and Path B (K > 1) the
  # parameters{} block of the inherited template differs structurally
  # (a_raw/c_b_raw flat-packed, theta_ref vector valued, K-slot
  # variants), so we delegate the init to cmdstanr's default
  # unconstrained-space random sampler (one per chain; reproducible
  # via the per-init seed offset). Multi-start coverage is preserved
  # because each iteration uses a distinct seed_offset.
  is_multi_or_K <- (isTRUE(stan_data$p_dim > 1L) ||
                    isTRUE(stan_data$K_slots > 1L))
  for (m in seq_len(M)) {
    init_m <- if (is_multi_or_K) NULL else {
      .gdpar_eb_make_random_init(stan_data, seed_offset = m,
                                 base_seed = seed)
    }
    # cmdstanr expects init as a list of lists (one per chain); we
    # treat each multi-start init as a single chain. jacobian = TRUE
    # is required for the downstream laplace() call to match the
    # unconstrained-scale convention.
    optim_args <- list(
      data      = stan_data,
      refresh   = 0L,
      algorithm = control$optim_algorithm,
      jacobian  = TRUE
    )
    if (!is.null(init_m)) {
      optim_args$init <- list(init_m)
    }
    if (!is.null(seed)) optim_args$seed <- as.integer(seed) + m
    opt_m <- tryCatch(
      do.call(model$optimize, optim_args),
      error = function(e) {
        if (verbose) {
          gdpar_inform(
            sprintf("Multi-start optimize() init %d failed: %s",
                    m, conditionMessage(e)),
            class = "gdpar_eb_message"
          )
        }
        NULL
      }
    )
    if (is.null(opt_m)) next
    lp_m <- tryCatch(opt_m$mle()["lp__"], error = function(e) NA_real_)
    history_lp[m] <- as.numeric(lp_m)
    if (is.na(best_idx) ||
        (isTRUE(is.finite(lp_m)) &&
         (is.na(history_lp[best_idx]) || lp_m > history_lp[best_idx]))) {
      best_idx <- m
      best_opt <- opt_m
      best_init <- init_m
    }
  }
  if (is.null(best_opt)) {
    gdpar_abort(
      paste(
        "All multi_start_M = ", M, " Laplace optimize() inits failed.",
        "Recommend using gdpar() (FB) on this model.", sep = ""
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(history_lp = history_lp)
    )
  }
  laplace_args <- list(
    data       = stan_data,
    mode       = best_opt,
    jacobian   = TRUE,
    draws      = control$laplace_draws,
    refresh    = 0L
  )
  if (!is.null(seed)) laplace_args$seed <- as.integer(seed) + 1000L
  lap <- tryCatch(
    do.call(model$laplace, laplace_args),
    error = function(e) NULL
  )
  if (is.null(lap)) {
    gdpar_abort(
      paste(
        "cmdstanr::laplace() failed on the best init returned by",
        "multi_start optimize(); Hessian likely singular or",
        "non-positive-definite at the candidate MAP."
      ),
      class = "gdpar_eb_numerical_error",
      data = list(history_lp = history_lp, best_idx = best_idx)
    )
  }
  draws <- lap$draws()
  # Path-aware theta_ref variable extraction (Sub-phase 8.6.C, D35):
  # K = 1 + p = 1 (8.6.B) -> "theta_ref" (J = 1) or "theta_ref[1..J]".
  # Path A K = 1 + p > 1   -> "theta_ref[j,k]" for j in 1..J, k in 1..p.
  # Path B K > 1 + p = 1   -> "theta_ref_k[j,k]" for j in 1..J, k in 1..K.
  all_vars <- dimnames(draws)$variable %||% character(0L)
  is_path_B <- isTRUE(stan_data$K_slots > 1L)
  is_path_A <- (!is_path_B) && isTRUE(stan_data$p_dim > 1L)
  if (is_path_B) {
    K_dim <- as.integer(stan_data$K_slots)
    n_anchors <- J * K_dim
    available <- character(0L)
    for (j in seq_len(J)) {
      for (k in seq_len(K_dim)) {
        available <- c(available, sprintf("theta_ref_k[%d,%d]", j, k))
      }
    }
    available <- intersect(available, all_vars)
    if (length(available) != n_anchors) {
      gdpar_abort(
        sprintf(
          "Path B EB Laplace draws missing theta_ref_k entries: expected %d found %d.",
          n_anchors, length(available)
        ),
        class = "gdpar_internal_error"
      )
    }
  } else if (is_path_A) {
    p_dim <- as.integer(stan_data$p_dim)
    n_anchors <- J * p_dim
    available <- character(0L)
    for (j in seq_len(J)) {
      for (k in seq_len(p_dim)) {
        available <- c(available, sprintf("theta_ref[%d,%d]", j, k))
      }
    }
    available <- intersect(available, all_vars)
    if (length(available) != n_anchors) {
      gdpar_abort(
        sprintf(
          "Path A EB Laplace draws missing theta_ref entries: expected %d found %d.",
          n_anchors, length(available)
        ),
        class = "gdpar_internal_error"
      )
    }
  } else {
    theta_vars <- paste0("theta_ref[", seq_len(J), "]")
    available <- intersect(theta_vars, all_vars)
    if (length(available) != J) {
      if (J == 1L) {
        single_var <- intersect("theta_ref", all_vars)
        if (length(single_var) == 1L) available <- single_var
      }
    }
    if (length(available) == 0L) {
      gdpar_abort(
        "theta_ref variable not found in Laplace draws output.",
        class = "gdpar_internal_error"
      )
    }
  }
  theta_draws <- posterior::subset_draws(draws, variable = available)
  theta_mat <- posterior::as_draws_matrix(theta_draws)
  theta_hat <- colMeans(theta_mat)
  theta_cov <- if (ncol(theta_mat) > 1L) stats::cov(theta_mat)
               else matrix(stats::var(theta_mat[, 1L]), 1L, 1L)
  # O5*-EBFB anti-fragility (Sub-bloque 9.3.b, B9.2): delegate the
  # adaptive Levenberg-Marquardt ridge to the canonical helper. The
  # helper triggers on non-PD eigenstructure or |det(cov)| < epsilon_lm,
  # then grows lambda geometrically until kappa_post <= kappa_threshold
  # or ridge_max_iter is reached. The post-ridge kappa is what gets
  # gated against kappa_threshold below — preserving the canonical
  # §2.8 semantics that the threshold applies to the (possibly ridged)
  # final covariance.
  lm_out <- .gdpar_eb_lm_perturb(theta_cov, control)
  theta_cov <- lm_out$cov_perturbed
  ridge_used <- lm_out$lambda_used
  kappa_val <- lm_out$kappa_post
  eigs <- tryCatch(
    eigen(theta_cov, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (isTRUE(kappa_val > control$kappa_threshold) ||
      identical(lm_out$status, "exhausted")) {
    gdpar_abort(
      sprintf(
        paste(
          "Marginal Hessian condition number kappa = %.3g exceeds",
          "threshold %.3g (laplace_control$kappa_threshold) after",
          "adaptive Levenberg-Marquardt ridge (status = '%s',",
          "n_iter = %d, lambda = %.3g). Multivariate Laplace is",
          "unreliable for this configuration; consider gdpar() (FB)",
          "instead. Smallest eigenvalue: %.3g."
        ),
        kappa_val, control$kappa_threshold, lm_out$status,
        lm_out$n_iter, ridge_used, min(eigs)
      ),
      class = "gdpar_eb_numerical_error",
      data = list(kappa = kappa_val, eigenvalues = eigs,
                  history_lp = history_lp,
                  lm_status = lm_out$status,
                  lm_n_iter = lm_out$n_iter,
                  lm_lambda = ridge_used)
    )
  }
  theta_se <- sqrt(pmax(diag(theta_cov), 0))
  history_clean <- history_lp[is.finite(history_lp)]
  dispersion <- if (length(history_clean) > 1L) {
    sd(history_clean) / max(abs(mean(history_clean)), 1)
  } else {
    NA_real_
  }
  if (isTRUE(dispersion > 0.05) && verbose) {
    gdpar_warn(
      sprintf(
        paste(
          "Multi-start dispersion across %d Laplace inits is %.3f",
          "(threshold 0.05): possible multimodality of the marginal",
          "likelihood (v07b Section 9.5 O5*-EBFB). Consider gdpar() (FB)."
        ),
        sum(is.finite(history_lp)), dispersion
      ),
      class = "gdpar_diagnostic_warning",
      data = list(dispersion = dispersion, history_lp = history_lp)
    )
  }
  list(
    theta_ref_hat = unname(theta_hat),
    theta_ref_se  = unname(theta_se),
    theta_ref_cov = theta_cov,
    diagnostics = list(
      kappa                       = kappa_val,
      lm_perturbation             = ridge_used,
      lm_n_iter                   = lm_out$n_iter,
      lm_status                   = lm_out$status,
      kappa_post_ridge            = lm_out$kappa_post,
      multi_start_dispersion      = dispersion,
      marginal_log_lik_history    = history_lp,
      best_init_index             = best_idx
    )
  )
}

#' Internal: random init helper for multi-start Laplace
#'
#' Produces a small named list passable to
#' \code{cmdstanr::optimize(init = ...)}. Inits are sampled from
#' standard normals for unconstrained parameters; \code{seed_offset}
#' lets the caller request reproducible inits.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_make_random_init <- function(stan_data, seed_offset = 1L,
                                        base_seed = NULL) {
  rng_seed <- if (is.null(base_seed)) NULL
              else as.integer(base_seed) + seed_offset
  if (!is.null(rng_seed)) {
    old_state <- if (exists(".Random.seed", envir = .GlobalEnv))
                   get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(rng_seed)
    on.exit({
      if (!is.null(old_state)) {
        assign(".Random.seed", old_state, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  J <- stan_data$J_groups
  init <- list()
  init$theta_ref <- stats::rnorm(J, mean = 0, sd = 0.5)
  if (isTRUE(stan_data$use_groups == 1L)) {
    init$mu_theta_ref    <- as.array(stats::rnorm(1L, sd = 0.5))
    init$sigma_theta_ref <- as.array(abs(stats::rnorm(1L)) + 0.1)
  }
  if (isTRUE(stan_data$use_a == 1L) && stan_data$J_a > 0L) {
    init$sigma_a <- as.array(abs(stats::rnorm(1L)) + 0.1)
    init$a_raw   <- stats::rnorm(stan_data$J_a, sd = 0.5)
  }
  if (isTRUE(stan_data$use_b == 1L) && stan_data$J_b > 0L) {
    init$sigma_b <- as.array(abs(stats::rnorm(1L)) + 0.1)
    init$c_b_raw <- stats::rnorm(stan_data$J_b, sd = 0.5)
  }
  if (isTRUE(stan_data$use_W == 1L) && stan_data$dim_W > 0L) {
    init$sigma_W <- as.array(abs(stats::rnorm(1L)) + 0.1)
    init$W_raw   <- matrix(stats::rnorm(stan_data$dim_W * stan_data$d, sd = 0.5),
                            nrow = stan_data$dim_W,
                            ncol = stan_data$d)
  }
  if (isTRUE(stan_data$use_dispersion_y == 1L)) {
    init$sigma_y <- as.array(abs(stats::rnorm(1L)) + 0.1)
  }
  if (isTRUE(stan_data$use_dispersion_phi == 1L)) {
    init$phi <- as.array(abs(stats::rnorm(1L)) + 1)
  }
  init
}

#' Internal: apply the Proposition 7B coverage-discrepancy correction
#'
#' Sub-phase 8.6.B: scalar form under p = 1, K = 1
#'   C_{g,alpha} = kappa(alpha) * (g'(xi*))^2 * (J^xi)^2 /
#'                 I_{theta theta}^{marg}
#' of v07 Section 6 (and v07b Section 5.2's reduction to scalar).
#'
#' Sub-phase 8.6.C (decision D34, Path A): when \code{p > 1} the
#' correction is matrix-valued (Proposition 7B*, v07b Section 5.1).
#' Delegates to \code{.gdpar_eb_correction_matrix()} which returns a
#' \eqn{p \times p} matrix \eqn{C^*_{g, \alpha} = \kappa(\alpha) \cdot
#' \Sigma^{marg}_{\theta_{ref}}} and reduces algebraically to the
#' scalar of 8.6.B at \eqn{p = 1}.
#'
#' For the default identity functional g(xi) = xi and a smooth model,
#' the framework reports the scaling object only; the actual inflation
#' of the credible intervals is applied downstream at the S3
#' \code{summary()} / \code{predict()} methods so the raw conditional
#' draws remain the canonical posterior representation.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_apply_correction <- function(eb_correction, laplace_result,
                                       stan_data, p = 1L, verbose) {
  if (!isTRUE(eb_correction)) {
    if (verbose) {
      gdpar_warn(
        paste(
          "eb_correction = FALSE: credible intervals will use nominal",
          "coverage and may under-cover by O(n^-1) per Proposition 7B",
          "of v07 Section 6."
        ),
        class = "gdpar_diagnostic_warning"
      )
    }
    return(list(applied = FALSE, constant = NA_real_))
  }
  p <- as.integer(p)
  if (p > 1L) {
    return(.gdpar_eb_correction_matrix(eb_correction = eb_correction,
                                       laplace_result = laplace_result,
                                       stan_data = stan_data,
                                       p = p,
                                       verbose = verbose))
  }
  kappa_alpha_95 <- 1.92
  marg_var <- tryCatch(
    laplace_result$theta_ref_cov[1L, 1L],
    error = function(e) NA_real_
  )
  if (!isTRUE(is.finite(marg_var)) || marg_var <= 0) {
    return(list(applied = FALSE, constant = NA_real_))
  }
  C_scalar <- kappa_alpha_95 * marg_var
  list(applied = TRUE, constant = as.double(C_scalar))
}

#' Internal: Proposition 7B* matrix-valued coverage-discrepancy correction
#'
#' Sub-phase 8.6.C (decision D34, Path A). Implements the matrix form
#' of v07b Section 5.1 Proposition 7B*:
#'   C^*_{g, alpha} = kappa(alpha) * J^xi^T * Sigma^marg_{theta_ref} * J^xi
#' For the default identity functional g(xi) = xi the Jacobian
#' \eqn{J^\xi} reduces to the identity and the correction equals
#' \eqn{\kappa(\alpha) \cdot \Sigma^{marg}_{\theta_{ref}}}, which is
#' the algebraic generalization of the scalar 8.6.B form: at
#' \eqn{p = 1} the \eqn{1 \times 1} matrix collapses to the scalar
#' \eqn{\kappa(\alpha) \cdot \mathrm{Var}^{marg}(\theta_{ref})}.
#'
#' Aborts the matrix-valued path silently to (applied = FALSE,
#' constant = NA matrix) when the marginal covariance is not finite or
#' not positive semi-definite; downstream S3 methods then fall back
#' to nominal credible intervals and a diagnostic warning is issued.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_correction_matrix <- function(eb_correction, laplace_result,
                                         stan_data, p = 1L, verbose) {
  p <- as.integer(p)
  if (!isTRUE(eb_correction)) {
    return(list(applied = FALSE, constant = matrix(NA_real_, p, p)))
  }
  kappa_alpha_95 <- 1.92
  cov_mat <- tryCatch(
    laplace_result$theta_ref_cov,
    error = function(e) NULL
  )
  if (is.null(cov_mat) || !is.matrix(cov_mat) ||
      nrow(cov_mat) != ncol(cov_mat)) {
    return(list(applied = FALSE, constant = matrix(NA_real_, p, p)))
  }
  if (any(!is.finite(cov_mat))) {
    return(list(applied = FALSE, constant = matrix(NA_real_, p, p)))
  }
  eigs <- tryCatch(eigen(cov_mat, symmetric = TRUE,
                          only.values = TRUE)$values,
                    error = function(e) NA_real_)
  if (any(!is.finite(eigs)) || any(eigs < -1e-10)) {
    return(list(applied = FALSE,
                constant = matrix(NA_real_, nrow(cov_mat), ncol(cov_mat))))
  }
  C_mat <- kappa_alpha_95 * cov_mat
  list(applied = TRUE, constant = C_mat)
}

#' Internal: resolve K-input patterns into amm_list_canonical
#'
#' Sub-phase 8.6.C (decision D34): mirrors the K-input dispatch of
#' \code{gdpar()} (gdpar.R section preceding \code{.gdpar_K}). Builds
#' the canonical \code{amm_list_canonical} from one of three input
#' shapes (\code{gdpar_formula_set}; named list of \code{amm_spec};
#' legacy formula with \code{a()}/\code{b()}/\code{W()} wrappers in
#' the RHS), promotes the family scope per slot, and reports the
#' resolved \code{K} (length of the canonical list) plus the metadata
#' needed by \code{.gdpar_eb_run_K()} when K > 1.
#'
#' The logic is duplicated from \code{gdpar()} rather than refactored
#' into a shared helper to preserve the bit-exact behaviour of every
#' FB golden under Sub-phase 8.6.C; unification into a shared helper
#' is registered as a follow-on debt.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_resolve_K_inputs <- function(formula, amm, W, family,
                                       formula_set_input,
                                       amm_list_input,
                                       classic_with_amm_calls,
                                       family_is_named_list) {
  if (formula_set_input) {
    if (!.gdpar_is_default_amm_spec(amm)) {
      gdpar_abort(
        paste(
          "When 'formula' is a gdpar_formula_set, the 'amm' argument",
          "must remain at its default amm_spec(): the formula set is",
          "the canonical declaration of every AMM component per slot."
        ),
        class = "gdpar_input_error"
      )
    }
    amm_list_canonical <- .gdpar_formula_set_to_amm_spec_list(formula, W)
    outcome_name <- formula$outcome
    formula_env <- formula$env
  } else if (amm_list_input) {
    amm_list_canonical <- amm
    for (slot in names(amm_list_canonical)) {
      if (!nzchar(slot)) {
        gdpar_abort(
          paste(
            "Every entry of the named list passed to 'amm' must carry",
            "a non-empty slot name (the canonical parameter name)."
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
          "When 'amm' is a named list of amm_spec, 'formula' must be",
          "a two-sided formula 'y ~ ...'."
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
    first_eligible <- if (family_is_named_list) {
      family[[1L]]$param_specs[[1L]]$name
    } else {
      family$param_specs[[1L]]$name
    }
    args_for_fs <- list(formula)
    names(args_for_fs) <- first_eligible
    fs <- do.call(gdpar_formula_set, args_for_fs)
    amm_list_canonical <- .gdpar_formula_set_to_amm_spec_list(fs, W)
    outcome_name <- fs$outcome
    formula_env <- fs$env
  }
  K <- length(amm_list_canonical)
  if (K > 1L) {
    if (family_is_named_list) {
      het <- .gdpar_resolve_heterogeneous_family_K(
        family, names(amm_list_canonical)
      )
      family_promoted     <- het$location_family
      family_id_k_vector  <- het$family_id_k_vector
    } else {
      family_promoted     <- .gdpar_promote_scope_per_observation(
        family, names(amm_list_canonical)
      )
      family_id_k_vector  <- NULL
    }
  } else {
    if (family_is_named_list) {
      gdpar_abort(
        paste(
          "A named list was supplied as 'family' (heterogeneous K-slot",
          "pattern) but the resolved K equals 1. The heterogeneous",
          "path requires K >= 2 K-individual slots; pass a single",
          "gdpar_family object for the K = 1 regime."
        ),
        class = "gdpar_input_error",
        data = list(K = K)
      )
    }
    k_name <- names(amm_list_canonical)[[1L]]
    family_promoted    <- .gdpar_promote_scope_per_observation(
      family, k_name
    )
    family_id_k_vector <- NULL
  }
  list(
    amm_list_canonical = amm_list_canonical,
    K                  = K,
    outcome_name       = outcome_name,
    formula_env        = formula_env,
    family_promoted    = family_promoted,
    family_id_k_vector = family_id_k_vector
  )
}

#' Internal: Path B EB orchestrator (K >= 2 distributional regression slots)
#'
#' Sub-phase 8.6.C (decision D34, Path B): companion of \code{.gdpar_K}
#' for the EB workflow. Reproduces the per-slot AMM design construction
#' and identifiability checks of \code{.gdpar_K}, then runs the
#' two-step EB pipeline: (i) Laplace approximation on
#' \code{amm_eb_marginal_K.stan} to obtain the per-slot EB anchors
#' \eqn{\widehat\theta_{ref,k}^{EB}}; (iii) HMC sampling on
#' \code{amm_eb_conditional_K.stan} with the plug-in anchors loaded as
#' data. Returns an object of class \code{gdpar_eb_fit} with the same
#' shape as the K = 1 path so the existing S3 methods accept either
#' regime uniformly.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_run_K <- function(amm_list_canonical, family, data, prior,
                            anchor, outcome_name, formula_env,
                            family_id_k_vector,
                            skip_id_check, chains, iter_warmup,
                            iter_sampling, adapt_delta,
                            max_treedepth, refresh, verbose, seed,
                            group, parametrization,
                            id_check_rigor, eb_correction,
                            laplace_control, call, ...) {
  slot_names <- names(amm_list_canonical)
  K <- length(amm_list_canonical)
  .gdpar_eb_check_stan_id_for_path(family, K, 1L)
  if (!outcome_name %in% colnames(data)) {
    gdpar_abort(
      sprintf("Outcome variable '%s' not found in data.", outcome_name),
      class = "gdpar_input_error"
    )
  }
  y <- data[[outcome_name]]
  if (is.matrix(y) || (is.array(y) && length(dim(y)) > 1L)) {
    gdpar_abort(
      sprintf(
        paste(
          "When K > 1 with p = 1, outcome '%s' must be a length-n",
          "vector (one univariate response shared across the K slots)."
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
        "Outcome '%s' contains %d non-finite value(s) (NA, NaN or Inf).",
        outcome_name, sum(y_bad)
      ),
      class = "gdpar_input_error"
    )
  }
  union_vars <- character(0L)
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]
    if (!is.null(a_k$a)) union_vars <- c(union_vars, all.vars(a_k$a))
    if (!is.null(a_k$b)) union_vars <- c(union_vars, all.vars(a_k$b))
  }
  union_vars <- unique(union_vars)
  rhs_str <- if (length(union_vars) == 0L) "1" else
             paste(union_vars, collapse = " + ")
  formula_full <- stats::as.formula(
    paste(outcome_name, "~", rhs_str),
    env = formula_env
  )
  rhs <- formula_full[c(1L, 3L)]
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
        if (!is.null(a_k$a)) all.vars(a_k$a) else character(0L),
        if (!is.null(a_k$b)) all.vars(a_k$b) else character(0L)
      ))
      slot_rhs_str <- if (length(slot_vars) == 0L) "1" else
                       paste(slot_vars, collapse = " + ")
      slot_rhs <- stats::as.formula(
        paste("~", slot_rhs_str), env = formula_env
      )
      rep_k <- gdpar_check_identifiability(
        amm = a_k, data = data,
        theta_ref_init = diag_theta,
        formula_rhs = slot_rhs
      )
      if (!isTRUE(rep_k$passed)) {
        gdpar_abort(
          sprintf(
            paste(
              "Basis-restricted identifiability check failed for slot",
              "'%s' at the diagnostic test point. Inspect the returned",
              "report or pass skip_id_check = TRUE to override at your",
              "own risk."
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
      gdpar_abort(
        paste(
          "K-level identifiability check failed. Inspect the returned",
          "report or pass skip_id_check = TRUE to override at your",
          "own risk."
        ),
        class = "gdpar_identifiability_error",
        data = list(report = K_level)
      )
    }
    attr(per_slot_reports, "K_level") <- K_level
    per_slot_reports
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
  stan_data$K_slots <- as.integer(K)
  stan_data$p_dim   <- 1L
  cp_a <- identical(parametrization, "cp")
  cp_W <- identical(parametrization, "cp")
  parametrization_resolved <- list(
    cp_a       = isTRUE(cp_a),
    cp_W       = isTRUE(cp_W),
    cp_a_per_K = NULL,
    meta = list(
      mode = "eb_K_path_B",
      note = paste(
        "Sub-phase 8.6.C Path B ships with uniform CP/NCP across the",
        "K slots; per-slot preflight (cp_a_per_K) is queued."
      ),
      requested = list(parametrization = parametrization)
    )
  )
  marginal_src <- .gdpar_eb_generate_stan_marginal(
    prior,
    cp_a   = parametrization_resolved$cp_a,
    cp_W   = parametrization_resolved$cp_W,
    K      = K,
    p      = 1L,
    family = family
  )
  marginal_path <- write_stan_to_tempfile(marginal_src)
  require_suggested("cmdstanr",
                    "fit Path 1 Empirical-Bayes models via cmdstanr::laplace()")
  marginal_model <- cmdstanr::cmdstan_model(marginal_path)
  laplace_result <- .gdpar_eb_maximize_marginal(
    model     = marginal_model,
    stan_data = stan_data,
    control   = laplace_control,
    seed      = seed,
    verbose   = verbose
  )
  conditional_src <- .gdpar_eb_generate_stan_conditional(
    prior,
    cp_a   = parametrization_resolved$cp_a,
    cp_W   = parametrization_resolved$cp_W,
    K      = K,
    p      = 1L,
    family = family
  )
  conditional_path  <- write_stan_to_tempfile(conditional_src)
  conditional_model <- cmdstanr::cmdstan_model(conditional_path)
  stan_data_cond <- stan_data
  # theta_ref_k_data: per-group array[J_groups] vector[K]. The Laplace
  # helper returns a flat numeric vector of length J_groups * K (in
  # row-major order over groups). For the default J_groups == 1 the
  # vector has length K and is reshaped to a 1-by-K matrix; the Stan
  # array[J_groups] vector[K] consumer accepts both R matrices and
  # nested lists of vectors via cmdstanr's automatic packing.
  theta_hat <- as.numeric(laplace_result$theta_ref_hat)
  J_groups_loc <- as.integer(stan_data$J_groups)
  if (length(theta_hat) == J_groups_loc * K) {
    stan_data_cond$theta_ref_k_data <- matrix(theta_hat,
                                              nrow = J_groups_loc,
                                              ncol = K, byrow = FALSE)
  } else {
    stan_data_cond$theta_ref_k_data <- matrix(theta_hat,
                                              nrow = J_groups_loc,
                                              ncol = K, byrow = FALSE)
  }
  sample_args <- list(
    data            = stan_data_cond,
    chains          = as.integer(chains),
    iter_warmup     = as.integer(iter_warmup),
    iter_sampling   = as.integer(iter_sampling),
    adapt_delta     = adapt_delta,
    max_treedepth   = as.integer(max_treedepth),
    refresh         = as.integer(refresh),
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) sample_args$seed <- as.integer(seed)
  extra <- list(...)
  for (nm in names(extra)) sample_args[[nm]] <- extra[[nm]]
  fit_cond <- do.call(conditional_model$sample, sample_args)
  diagnostics <- compute_diagnostics(fit_cond, verbose = verbose)
  correction <- .gdpar_eb_apply_correction(
    eb_correction  = eb_correction,
    laplace_result = laplace_result,
    stan_data      = stan_data,
    p              = 1L,
    verbose        = verbose
  )
  obj <- list(
    theta_ref_hat            = laplace_result$theta_ref_hat,
    theta_ref_se             = laplace_result$theta_ref_se,
    conditional_fit          = fit_cond,
    amm_list_canonical       = amm_list_canonical,
    family                   = family,
    prior                    = prior,
    design_K                 = design_K,
    anchor                   = anchor_value,
    stan_data                = stan_data,
    identifiability_report   = id_report,
    diagnostics              = diagnostics,
    diagnostics_numerical    = laplace_result$diagnostics,
    parametrization          = parametrization_resolved,
    group_info               = group_info,
    correction_applied       = correction$applied,
    eb_correction_constant   = correction$constant,
    call                     = call,
    path                     = "eb",
    K                        = K,
    slot_names               = slot_names
  )
  class(obj) <- c("gdpar_eb_fit", "list")
  obj
}

#' Internal: Proposition 7B* tensor-valued coverage-discrepancy correction
#'
#' Sub-phase 8.6.D (decision D37 = (i) canonized Session 13a 2026-05-25
#' + decision D38'' = (h) canonized Session 13b 2026-05-25). Implements
#' the tensor form of v07b Section 5.1 Proposition 7B* under the full
#' regime K > 1 AND p > 1:
#'   C^{tensor}_{g, alpha}[k, , ] = kappa(alpha) * Sigma^{marg}_{theta_ref, k}
#' where \eqn{\Sigma^{marg}_{\theta_{ref}, k} \in \mathbb{R}^{p \times p}}
#' is the marginal covariance matrix of the per-slot anchor
#' \eqn{\theta_{ref, k}} across its p coordinates, retaining the
#' cross-coordinate terms inside each slot. Total shape: 3D array
#' \code{[K, p, p]} with all cross terms inside each slot retained
#' (D37 = (i) "tensor full"). The block-diagonal structure across
#' slots is the canonical assumption of the prior factorization
#' \eqn{\pi_\xi = \prod_k \pi_{\xi_k}}; under coupling priors the
#' off-block terms are zero by construction of this estimator (the
#' Laplace marginal covariance is computed slot-by-slot in the Path C
#' orchestrator).
#'
#' At p = 1 the tensor degenerates to a \code{[K, 1, 1]} array whose
#' slot blocks recover the scalar Path B Proposition 7B correction
#' coordinate by coordinate. At K = 1 the tensor degenerates to a
#' \code{[1, p, p]} array whose single block is the Path A matrix
#' correction.
#'
#' Aborts the tensor-valued path silently to (applied = FALSE,
#' constant = NA tensor) when any per-slot covariance is not finite,
#' not positive semi-definite, or unavailable; downstream S3 methods
#' then fall back to nominal credible intervals and a diagnostic
#' warning is issued.
#'
#' @param eb_correction Logical scalar.
#' @param laplace_result_per_slot Named list with K entries; each entry
#'   carries the slot-specific Laplace covariance \eqn{p \times p}
#'   matrix as \code{theta_ref_cov_k}.
#' @param K Integer scalar; number of slots.
#' @param p Integer scalar; outcome dimension.
#' @param verbose Logical scalar.
#'
#' @return A list with \code{applied} (logical), \code{constant} (3D
#'   array \code{[K, p, p]} with the per-slot correction blocks; NA
#'   filled when applied = FALSE), and \code{slot_dispositions} (named
#'   character with the per-slot disposition: "ok", "non_finite",
#'   "non_psd", "missing").
#'
#' @keywords internal
#' @noRd
.gdpar_eb_correction_tensor <- function(eb_correction,
                                         laplace_result_per_slot,
                                         K = 2L, p = 1L, verbose = TRUE) {
  K <- as.integer(K)
  p <- as.integer(p)
  empty_tensor <- function() {
    arr <- array(NA_real_, dim = c(K, p, p))
    arr
  }
  if (!isTRUE(eb_correction)) {
    return(list(applied            = FALSE,
                constant           = empty_tensor(),
                slot_dispositions  = setNames(rep("disabled", K),
                                              seq_len(K))))
  }
  kappa_alpha_95 <- 1.92
  tensor <- empty_tensor()
  dispositions <- character(K)
  any_failed <- FALSE
  for (k in seq_len(K)) {
    slot_lap <- laplace_result_per_slot[[k]]
    cov_k <- tryCatch(slot_lap$theta_ref_cov_k, error = function(e) NULL)
    if (is.null(cov_k)) {
      dispositions[k] <- "missing"
      any_failed <- TRUE
      next
    }
    if (!is.matrix(cov_k) || nrow(cov_k) != p || ncol(cov_k) != p) {
      dispositions[k] <- "missing"
      any_failed <- TRUE
      next
    }
    if (any(!is.finite(cov_k))) {
      dispositions[k] <- "non_finite"
      any_failed <- TRUE
      next
    }
    eigs_k <- tryCatch(
      eigen(cov_k, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    if (any(!is.finite(eigs_k)) || any(eigs_k < -1e-10)) {
      dispositions[k] <- "non_psd"
      any_failed <- TRUE
      next
    }
    tensor[k, , ] <- kappa_alpha_95 * cov_k
    dispositions[k] <- "ok"
  }
  if (any_failed) {
    if (verbose) {
      gdpar_warn(
        sprintf(
          paste(
            "Path C correction tensor: %d of %d slots failed",
            "(%s); downstream S3 methods fall back to nominal coverage."
          ),
          sum(dispositions != "ok"), K,
          paste(setdiff(unique(dispositions), "ok"), collapse = ", ")
        ),
        class = "gdpar_diagnostic_warning"
      )
    }
    return(list(applied = FALSE, constant = empty_tensor(),
                slot_dispositions = dispositions))
  }
  list(applied = TRUE, constant = tensor,
       slot_dispositions = dispositions)
}

#' Internal: build the design matrices for a Path C K x p specification
#'
#' Sub-phase 8.6.D (decision D38'' = (h) + D41 canonized Session 13b
#' 2026-05-25 + Session 13c 2026-05-26 implementation). Path C consumes
#' a named list of \code{amm_spec} objects of length K where every slot
#' carries \code{p > 1}; this helper iterates over the K slots and per
#' slot delegates to \code{.build_amm_design_multi()} to produce the
#' per-coord ragged design matrices. The result is the input that
#' \code{.assemble_stan_data_KxP()} consumes to populate the data
#' fields of \code{amm_eb_marginal_KxP.stan}.
#'
#' @param amm_list_canonical Named list of K amm_spec objects with
#'   \code{p > 1} per slot (homogeneous p enforced upstream by the
#'   orchestrator).
#' @param data Data frame containing the variables referenced by the
#'   per-slot amm specifications.
#' @param formula_rhs Two-sided formula identifying the covariate
#'   columns of \code{data} used as the linear factor x; passed through
#'   to \code{.build_amm_design_multi()} per slot.
#'
#' @return Named list with K, p, slot_names, and three K-length lists
#'   of per-slot multivariate design objects: \code{design_per_slot}
#'   (each entry is the list returned by .build_amm_design_multi()).
#' @keywords internal
#' @noRd
.build_amm_design_KxP <- function(amm_list_canonical, data, formula_rhs) {
  if (!is.list(amm_list_canonical) || length(amm_list_canonical) < 2L) {
    gdpar_abort(
      "Internal error: .build_amm_design_KxP expected a named list of length K >= 2.",
      class = "gdpar_internal_error"
    )
  }
  slot_names <- names(amm_list_canonical)
  if (is.null(slot_names) || any(!nzchar(slot_names))) {
    gdpar_abort(
      "Internal error: .build_amm_design_KxP expected every slot to carry a non-empty name.",
      class = "gdpar_internal_error"
    )
  }
  K <- length(amm_list_canonical)
  p_per_slot <- vapply(amm_list_canonical,
                        function(a) as.integer(a$p %||% 1L),
                        integer(1L))
  if (any(p_per_slot < 2L) || any(p_per_slot != p_per_slot[1L])) {
    gdpar_abort(
      sprintf(
        paste(
          "Internal error: .build_amm_design_KxP requires homogeneous p >= 2",
          "across slots; received p_per_slot = %s."
        ),
        paste(p_per_slot, collapse = ", ")
      ),
      class = "gdpar_internal_error"
    )
  }
  p <- p_per_slot[1L]
  assert_data_frame(data, "data")
  design_per_slot <- vector("list", K)
  names(design_per_slot) <- slot_names
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]
    assert_inherits(a_k, "amm_spec", paste0("amm_list_canonical[['",
                                             slot_names[k], "']]"))
    design_per_slot[[k]] <- .build_amm_design_multi(a_k, data, formula_rhs)
  }
  list(
    K               = K,
    p               = p,
    slot_names      = slot_names,
    design_per_slot = design_per_slot
  )
}

#' Internal: assemble the data block passed to the Path C Stan template
#'
#' Sub-phase 8.6.D (decisions D36 = alpha + D38'' = (h) + D39 + D40' +
#' D41 + D43 = (a)). Companion of \code{.assemble_stan_data_K} for the
#' Path C K x p regime. Consumes the design returned by
#' \code{.build_amm_design_KxP()} together with the (already canonized)
#' family and outcome matrix and produces the named list expected by
#' \code{amm_eb_marginal_KxP.stan} and \code{amm_eb_conditional_KxP.stan}.
#'
#' Restrictions enforced (per D39 + D40' + D41 of the parcial):
#' \itemize{
#'   \item \code{use_W} hardcoded to 0 (D39): the modulating W is
#'     deferred to Block 9.x. The Path C dispatcher upstream rejects
#'     input that declares \code{W != NULL} on any slot.
#'   \item \code{family_id_k_vector} homogeneous per slot, with
#'     \code{stan_id in {1, 3}} (D40'). The location-slot stan_id is
#'     replicated K times because the K=2 Path B parametrization
#'     hardcoded in the template branches mu / dispersion by slot
#'     index, not by per-slot stan_id.
#'   \item \code{Z_a_kp[k, j]} and \code{Z_b_kp[k, j]} are per-slot
#'     per-coord ragged (D41), padded to the global \code{J_a_max} /
#'     \code{J_b_max} per the Stan template contract.
#' }
#'
#' @param design_KxP List returned by \code{.build_amm_design_KxP()}.
#' @param family The (already promoted) gdpar_family with
#'   \code{stan_id in {1, 3}}.
#' @param amm_list_canonical Named list of K amm_spec objects with
#'   p >= 2 per slot. Used to extract per-slot use_a / use_b flags.
#' @param y_matrix Numeric or integer matrix of outcomes (n by p).
#' @param theta_anchor_kp Numeric matrix K by p with the per-slot
#'   per-coord anchors on the linear-predictor scale.
#' @param group_id Optional integer vector of length n.
#' @return Named list with the entries expected by
#'   \code{amm_eb_marginal_KxP.stan}.
#' @keywords internal
#' @noRd
.assemble_stan_data_KxP <- function(design_KxP, family, amm_list_canonical,
                                     y_matrix, theta_anchor_kp,
                                     group_id = NULL,
                                     path = c("EB", "FB"),
                                     cp_W = FALSE) {
  # Sub-bloque 9.3.d colateral D69 (Bloque 9, Sesion B9.6, 2026-05-27)
  # under decision M.iv unified canonized at session open: single
  # assembler with dispatch on `path` in {"EB", "FB"}. The EB branch
  # preserves the Sub-fase 8.6.D first iteration restrictions
  # (use_W = 0 per D39 + stan_id in {1, 3} per D40') byte-identical so
  # legacy EB callers keep their bit-exact data list. The FB branch
  # lifts both restrictions: W is enabled (globally shared, K * p *
  # W_per_kj_dim block of W_raw rows) and the Path B family set
  # {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13} is accepted, in line with
  # the canonical piece amm_canonical_pmulti_KxP.stan (B9.5 J.iv.A)
  # which expects use_W + the W type metadata as additional data
  # fields. Default path = "EB" preserves backward compatibility with
  # the .gdpar_eb_fit_KxP caller (eb.R ~ line 2937).
  path <- match.arg(path)
  if (!is.list(design_KxP) || is.null(design_KxP$design_per_slot)) {
    gdpar_abort(
      "Internal error: .assemble_stan_data_KxP expected a design list from .build_amm_design_KxP().",
      class = "gdpar_internal_error"
    )
  }
  K <- as.integer(design_KxP$K)
  p <- as.integer(design_KxP$p)
  if (K < 2L || p < 2L) {
    gdpar_abort(
      sprintf(
        "Internal error: .assemble_stan_data_KxP requires K >= 2 AND p >= 2; received K = %d, p = %d.",
        K, p
      ),
      class = "gdpar_internal_error"
    )
  }
  assert_inherits(family, "gdpar_family", "family")
  stan_id <- as.integer(family$stan_id)
  if (path == "EB") {
    if (!stan_id %in% c(1L, 3L)) {
      gdpar_abort(
        sprintf(
          paste(
            "Path C (Sub-phase 8.6.D first iteration) only supports stan_id",
            "in {1 (Gaussian K=2), 3 (Negative Binomial K=2)} per D40';",
            "received stan_id = %d. The remaining Path B set {5, 6, 7, 8,",
            "9, 10, 11, 12, 13} is deferred to a later iteration of 8.6.D",
            "with the explicit numerical caveat of opening Section 6.1."
          ),
          stan_id
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(stan_id = stan_id, family_name = family$name)
      )
    }
  } else {
    # FB path (B9.6 M.iv unified): the canonical piece
    # amm_canonical_pmulti_KxP.stan covers Path B family set
    # {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13} per the J.iv.A canonization
    # of B9.5; reject everything else.
    if (!stan_id %in% c(1L, 3L, 5L, 6L, 7L, 8L, 9L,
                          10L, 11L, 12L, 13L)) {
      gdpar_abort(
        sprintf(
          paste(
            "FB Path C KxP (Sub-bloque 9.3.d) supports stan_id in",
            "{1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13} per the J.iv.A piece",
            "canonized in B9.5; received stan_id = %d (family '%s')."
          ),
          stan_id, family$name
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(stan_id = stan_id, family_name = family$name)
      )
    }
  }
  if (!is.matrix(y_matrix)) {
    gdpar_abort(
      paste(
        "Internal error: .assemble_stan_data_KxP expected y_matrix to be",
        "an n by p numeric matrix."
      ),
      class = "gdpar_internal_error",
      data = list(received_class = class(y_matrix))
    )
  }
  if (ncol(y_matrix) != p) {
    gdpar_abort(
      sprintf(
        "Internal error: y_matrix has %d columns but Path C resolved p = %d.",
        ncol(y_matrix), p
      ),
      class = "gdpar_internal_error"
    )
  }
  n <- nrow(y_matrix)
  # Sub-bloque 9.3.d D69 fix (Bloque 9, Sesion B9.6, 2026-05-27 emergente
  # M.iv): for path = "FB" the full Path B family set extends beyond the
  # 8.6.D EB-side restriction {1 gaussian, 3 NB}. Continuous families
  # (gaussian / beta / gamma / lognormal_loc_scale / student_t / tweedie)
  # populate y_real; count families (NB / ZIP / ZINB / hurdle_poisson /
  # hurdle_neg_binomial_2) populate y_int. Without this extension y_real
  # is zeroed-out for path = "FB" with stan_id in {5, 6, 7, 8, 9}, and
  # beta_proportion / gamma / lognormal / student_t / tweedie lpdf
  # collapse to log(0) at sample-time, blocking initialization.
  needs_real <- stan_id %in% c(1L, 5L, 6L, 7L, 8L, 9L)
  needs_int  <- stan_id %in% c(3L, 10L, 11L, 12L, 13L)
  if (needs_int) {
    if (any(!is.finite(y_matrix)) ||
        any(y_matrix < 0) ||
        any(y_matrix != floor(y_matrix))) {
      gdpar_abort(
        paste(
          "Path C Negative Binomial K=2 requires every entry of the n by p",
          "outcome matrix to be a non-negative integer count.",
          "Rescale, round or trim before fitting."
        ),
        class = "gdpar_input_error",
        data = list(family_name = family$name,
                    y_range = range(y_matrix, na.rm = TRUE))
      )
    }
  }
  if (needs_real) {
    if (any(!is.finite(y_matrix))) {
      gdpar_abort(
        paste(
          "Path C Gaussian K=2 requires every entry of the n by p",
          "outcome matrix to be finite."
        ),
        class = "gdpar_input_error",
        data = list(y_range = range(y_matrix, na.rm = TRUE))
      )
    }
  }
  y_real <- if (needs_real) matrix(as.double(y_matrix), n, p)
            else matrix(0.0, n, p)
  y_int <- if (needs_int) matrix(as.integer(y_matrix), n, p)
           else matrix(0L, n, p)
  grp <- .resolve_group_id(group_id, n)
  use_a_k <- vapply(amm_list_canonical,
                     function(a) as.integer(!is.null(a$a) ||
                                              (!is.null(a$dims) &&
                                                any(vapply(a$dims,
                                                            function(d) !is.null(d$a),
                                                            logical(1L))))),
                     integer(1L))
  use_b_k <- vapply(amm_list_canonical,
                     function(a) as.integer(!is.null(a$b) ||
                                              (!is.null(a$dims) &&
                                                any(vapply(a$dims,
                                                            function(d) !is.null(d$b),
                                                            logical(1L))))),
                     integer(1L))
  any_W <- any(vapply(amm_list_canonical,
                       function(a) !is.null(a$W),
                       logical(1L)))
  if (path == "EB" && isTRUE(any_W)) {
    gdpar_abort(
      paste(
        "Path C (Sub-phase 8.6.D first iteration) has the modulating",
        "component W disabled per D39 on the EB path. At least one slot",
        "declares W != NULL; remove the W() wrapper or use path = 'FB'",
        "via the FB KxP fit harness (B9.6 D69 closure)."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  J_a_per_kp <- matrix(0L, nrow = K, ncol = p)
  J_b_per_kp <- matrix(0L, nrow = K, ncol = p)
  for (k in seq_len(K)) {
    d_k <- design_KxP$design_per_slot[[k]]
    for (j in seq_len(p)) {
      J_a_per_kp[k, j] <- ncol(d_k$Z_a_list[[j]])
      J_b_per_kp[k, j] <- ncol(d_k$Z_b_list[[j]])
    }
  }
  J_a_max <- max(J_a_per_kp, 0L)
  J_b_max <- max(J_b_per_kp, 0L)
  pad_to <- function(z, target_cols, n_rows) {
    if (target_cols == 0L) {
      return(matrix(0, nrow = n_rows, ncol = 0L))
    }
    if (ncol(z) >= target_cols) {
      return(z)
    }
    cbind(z, matrix(0, nrow = nrow(z), ncol = target_cols - ncol(z)))
  }
  # Z_a_kp / Z_b_kp packed as 4D arrays with shape [K, p, n, J_max] so
  # cmdstanr accepts them as array[K, p] matrix[n, J_max].
  Z_a_kp <- array(0.0, dim = c(K, p, n, J_a_max))
  Z_b_kp <- array(0.0, dim = c(K, p, n, J_b_max))
  for (k in seq_len(K)) {
    d_k <- design_KxP$design_per_slot[[k]]
    for (j in seq_len(p)) {
      z_a <- pad_to(d_k$Z_a_list[[j]], J_a_max, n)
      z_b <- pad_to(d_k$Z_b_list[[j]], J_b_max, n)
      if (J_a_max > 0L) Z_a_kp[k, j, , ] <- z_a
      if (J_b_max > 0L) Z_b_kp[k, j, , ] <- z_b
    }
  }
  # family_id_k_vector homogeneous per slot at the initial iteration
  # (D40'). The Stan template branches on family_id_k_vector[1] and
  # hardcodes the K=2 mu / dispersion semantics; the remaining slot
  # carries the same stan_id as a metadata convention.
  family_id_k_vector <- rep(as.integer(stan_id), K)
  inv_link_id_per_slot <- .gdpar_compute_inv_link_id_per_slot(
    family_id_k_vector, family
  )
  if (!is.matrix(theta_anchor_kp) ||
      nrow(theta_anchor_kp) != K ||
      ncol(theta_anchor_kp) != p) {
    gdpar_abort(
      sprintf(
        paste(
          "Internal error: theta_anchor_kp must be a K by p numeric matrix",
          "(K = %d, p = %d); received an object of dim %s."
        ),
        K, p,
        paste(if (is.null(dim(theta_anchor_kp))) length(theta_anchor_kp)
              else dim(theta_anchor_kp), collapse = " by ")
      ),
      class = "gdpar_internal_error"
    )
  }
  theta_anchor_kp_list <- lapply(seq_len(K),
                                  function(k) as.double(theta_anchor_kp[k, ]))
  if (path == "EB") {
    use_W_value <- 0L
    use_dispersion_y_k   <- as.integer(rep(0L, K))
    use_dispersion_phi_k <- as.integer(rep(0L, K))
  } else {
    # FB path (B9.6 M.iv): use_W = 1 if any slot declares W (W is
    # globally shared across slots and coords; the K * p * W_per_kj_dim
    # block of W_raw rows lives in the parameters{} block of
    # amm_canonical_pmulti_KxP.stan, indexed (k, j) via the
    # offset ((k-1)*p + (j-1)) * W_per_kj_dim + 1 per J.iv.A B9.5).
    use_W_value <- as.integer(any_W)
    # Per-slot population-scope dispersion flags: gaussian / NB define
    # sigma / phi as POPULATION-scope by default; B9.5 J.iv.A piece
    # accepts both K-individual (per-observation) and population-scope
    # dispersion via the sigma_y_pop_k / phi_pop_k arrays gated by
    # use_dispersion_y_k / use_dispersion_phi_k. For the B9.6 atomic
    # closure of D69, the harness uses K-individual dispersion via
    # explicit aux slots (mirror of the EB-side); the flags stay zero.
    # Future B9.7+ may lift this for the bench external 9.2 if needed.
    use_dispersion_y_k   <- as.integer(rep(0L, K))
    use_dispersion_phi_k <- as.integer(rep(0L, K))
  }
  base_list <- list(
    n                    = n,
    K                    = K,
    p                    = p,
    family_id_k_vector   = as.integer(family_id_k_vector),
    inv_link_id_per_slot = as.integer(inv_link_id_per_slot),
    use_a_k              = as.integer(use_a_k),
    use_b_k              = as.integer(use_b_k),
    use_W                = use_W_value,
    J_a_max              = as.integer(J_a_max),
    J_b_max              = as.integer(J_b_max),
    J_a_per_kp           = J_a_per_kp,
    J_b_per_kp           = J_b_per_kp,
    Z_a_kp               = Z_a_kp,
    Z_b_kp               = Z_b_kp,
    y_real               = y_real,
    y_int                = y_int,
    theta_anchor_kp      = theta_anchor_kp_list,
    use_dispersion_y_k   = use_dispersion_y_k,
    use_dispersion_phi_k = use_dispersion_phi_k,
    use_groups           = grp$use_groups,
    J_groups             = grp$J_groups,
    group_id             = grp$group_id,
    K_slots              = as.integer(K),
    p_dim                = as.integer(p)
  )
  if (path == "EB") {
    return(base_list)
  }
  # FB path: extend with W-block fields expected by
  # amm_canonical_pmulti_KxP.stan (X, dim_W, d, W_per_kj_dim,
  # W_type_id, W_n_knots_full, W_knots_full, W_degree). Same resolution
  # convention as .assemble_stan_data_multi: first slot that declares
  # W defines the basis metadata; all slots that declare W share that
  # canonical W object (enforced by the formula-set canonical
  # constructor upstream).
  first_W_amm <- NULL
  for (a in amm_list_canonical) {
    if (!is.null(a$W)) { first_W_amm <- a; break }
  }
  if (use_W_value == 1L && !is.null(first_W_amm)) {
    dim_W <- as.integer(first_W_amm$W$dim)
    # Globally shared W: total W_raw rows = K * p * W_per_kj_dim.
    if (dim_W <= 0L) {
      gdpar_abort(
        "Internal error: first slot's amm$W$dim must be positive when use_W = 1 on FB KxP.",
        class = "gdpar_internal_error",
        data = list(dim_W = dim_W)
      )
    }
    # `dim_W` per the W spec is interpreted as the per-(slot, coord)
    # basis dim in the canonical KxP piece; allow the natural
    # convention where amm$W$dim already gives W_per_kj_dim.
    W_per_kj_dim <- dim_W
    dim_W_total <- K * p * W_per_kj_dim
    X_full <- design_KxP$design_per_slot[[1L]]$X
    d <- if (is.null(X_full)) 0L else as.integer(ncol(X_full))
    X_mat <- if (d > 0L) X_full else matrix(0, nrow = n, ncol = 0L)
  } else {
    use_W_value <- 0L
    base_list$use_W <- 0L
    W_per_kj_dim <- 0L
    dim_W_total <- 0L
    d <- 0L
    X_mat <- matrix(0, nrow = n, ncol = 0L)
  }
  W_stan <- .gdpar_resolve_W_stan_data(
    W      = first_W_amm$W,
    use_W  = use_W_value,
    y_vec  = as.numeric(y_matrix),
    family = family
  )
  c(base_list, list(
    dim_W          = as.integer(dim_W_total),
    d              = as.integer(d),
    W_per_kj_dim   = as.integer(W_per_kj_dim),
    X              = X_mat,
    W_type_id      = W_stan$W_type_id,
    W_n_knots_full = W_stan$W_n_knots_full,
    W_knots_full   = W_stan$W_knots_full,
    W_degree       = W_stan$W_degree
  ))
}

#' Internal: random init helper for Path C multi-start Laplace
#'
#' Sub-phase 8.6.D (Session 13c 2026-05-26). The Path A / Path B
#' delegation to cmdstanr's default unconstrained-space random init
#' fails under Path C because the K x p geometry of the multi-parametric
#' likelihood explodes the joint when the random a_raw * sigma_a_k
#' product on slot 2 (dispersion) lands far from the origin (eta_kp[2]
#' >> 0 makes exp(eta_kp[2]) overflow to +Inf and breaks normal_lpdf).
#' This helper provides a conservative init: theta_ref_kp and a_raw at
#' 0 with small jitter; sigma_a_k / sigma_b_k at 0.1 + small jitter.
#' The init is path-specific and lives next to .gdpar_eb_make_random_init
#' to mirror the K = 1 + p = 1 helper while adapting to the K x p
#' parameter block of amm_eb_marginal_KxP.stan.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_make_random_init_KxP <- function(stan_data, seed_offset = 1L,
                                            base_seed = NULL) {
  rng_seed <- if (is.null(base_seed)) NULL
              else as.integer(base_seed) + seed_offset
  if (!is.null(rng_seed)) {
    old_state <- if (exists(".Random.seed", envir = .GlobalEnv))
                   get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(rng_seed)
    on.exit({
      if (!is.null(old_state)) {
        assign(".Random.seed", old_state, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  K <- as.integer(stan_data$K_slots %||% stan_data$K)
  p <- as.integer(stan_data$p_dim %||% stan_data$p)
  J <- as.integer(stan_data$J_groups %||% 1L)
  use_groups <- as.integer(stan_data$use_groups %||% 0L)
  any_use_a <- if (!is.null(stan_data$use_a_k))
                  as.integer(any(stan_data$use_a_k == 1L)) else 0L
  any_use_b <- if (!is.null(stan_data$use_b_k))
                  as.integer(any(stan_data$use_b_k == 1L)) else 0L
  init <- list()
  # theta_ref_kp: array[J, K] vector[p] -> 3D numeric array [J, K, p]
  # via cmdstanr's automatic packing convention (same as the conditional
  # template's theta_ref_kp_data plug-in).
  init$theta_ref_kp <- array(stats::rnorm(J * K * p, mean = 0, sd = 0.1),
                              dim = c(J, K, p))
  if (use_groups == 1L) {
    init$mu_theta_ref_kp <- array(
      stats::rnorm(1L * K * p, mean = 0, sd = 0.1),
      dim = c(1L, K, p)
    )
    init$sigma_theta_ref_kp <- array(
      abs(stats::rnorm(1L * K * p, mean = 0.5, sd = 0.05)),
      dim = c(1L, K, p)
    )
  }
  if (any_use_a == 1L) {
    # Option A (RG.6, D96): sigma_a_k is compacted to slots carrying
    # free a coefficients; mirror the n_sigma_a of the KxP template's
    # transformed data so the init vector matches the parameter length.
    # When every slot carries free a coefficients n_sigma_a == K, so the
    # rnorm draw count is unchanged and the init stays bit-identical.
    .jap <- matrix(as.integer(stan_data$J_a_per_kp), nrow = K, ncol = p)
    .slot_free_a <- (as.integer(stan_data$use_a_k) == 1L) &
                    (rowSums(.jap > 0L) > 0L)
    n_sigma_a <- sum(.slot_free_a)
    if (n_sigma_a > 0L) {
      init$sigma_a_k <- as.array(rep(0.1, n_sigma_a) +
                                   abs(stats::rnorm(n_sigma_a, sd = 0.02)))
    }
    total_J_a_free <- sum(stan_data$J_a_per_kp)
    if (total_J_a_free > 0L) {
      init$a_raw <- as.numeric(stats::rnorm(total_J_a_free, mean = 0,
                                              sd = 0.1))
    }
  }
  if (any_use_b == 1L) {
    init$sigma_b_k <- as.array(rep(0.1, K) +
                                 abs(stats::rnorm(K, sd = 0.02)))
    total_J_b_free <- sum(stan_data$J_b_per_kp)
    if (total_J_b_free > 0L) {
      init$c_b_kp_raw <- as.numeric(stats::rnorm(total_J_b_free, mean = 0,
                                                   sd = 0.1))
    }
  }
  init
}

#' Internal: Step (i) of the EB workflow under Path C
#'
#' Sub-phase 8.6.D (decision D43 = (a) canonized Session 13c
#' 2026-05-26). Path C maximizer mirroring \code{.gdpar_eb_maximize_marginal}
#' but specialized for the K x p regime: runs a single joint Laplace
#' approximation over the full tensor of \code{theta_ref_kp}
#' [J_groups x K x p] anchors (decision D43 = (a) joint Laplace) and
#' extracts the per-slot p by p covariance block from the marginal
#' draws (block-diagonal extraction post-hoc), packaging the result
#' into the per-slot list expected by
#' \code{.gdpar_eb_correction_tensor()}.
#'
#' The anti-fragility strategy of Charter Section 2.8 (multi_start,
#' ridge L-M, kappa check) is preserved: kappa is reported per slot,
#' and the slot-wise ridge L-M is applied lazily inside
#' \code{.gdpar_eb_correction_tensor()} via its own non-finite / non-PSD
#' branches.
#'
#' @return A list with components \code{theta_ref_kp_hat} (3D array
#'   [J_groups, K, p]), \code{laplace_result_per_slot} (named list of
#'   length K, each entry carries \code{theta_ref_cov_k} of shape p
#'   by p), \code{diagnostics} (kappa per slot, multi_start_dispersion,
#'   marginal_log_lik_history).
#' @keywords internal
#' @noRd
.gdpar_eb_maximize_marginal_KxP <- function(model, stan_data, control,
                                             seed, verbose) {
  M <- control$multi_start_M
  J <- stan_data$J_groups
  K <- stan_data$K_slots
  p <- stan_data$p_dim
  history_lp <- rep(NA_real_, M)
  best_idx <- NA_integer_
  best_opt <- NULL
  for (m in seq_len(M)) {
    init_m <- .gdpar_eb_make_random_init_KxP(stan_data, seed_offset = m,
                                              base_seed = seed)
    optim_args <- list(
      data      = stan_data,
      init      = list(init_m),
      refresh   = 0L,
      algorithm = control$optim_algorithm,
      jacobian  = TRUE
    )
    if (!is.null(seed)) optim_args$seed <- as.integer(seed) + m
    opt_m <- tryCatch(
      do.call(model$optimize, optim_args),
      error = function(e) {
        if (verbose) {
          gdpar_inform(
            sprintf("Path C multi-start optimize() init %d failed: %s",
                    m, conditionMessage(e)),
            class = "gdpar_eb_message"
          )
        }
        NULL
      }
    )
    if (is.null(opt_m)) next
    lp_m <- tryCatch(opt_m$mle()["lp__"], error = function(e) NA_real_)
    history_lp[m] <- as.numeric(lp_m)
    if (is.na(best_idx) ||
        (isTRUE(is.finite(lp_m)) &&
         (is.na(history_lp[best_idx]) || lp_m > history_lp[best_idx]))) {
      best_idx <- m
      best_opt <- opt_m
    }
  }
  if (is.null(best_opt)) {
    gdpar_abort(
      paste(
        "Path C: all multi_start_M = ", M, " Laplace optimize() inits failed.",
        " Recommend using gdpar() (FB) on this K x p configuration.",
        sep = ""
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(history_lp = history_lp)
    )
  }
  laplace_args <- list(
    data       = stan_data,
    mode       = best_opt,
    jacobian   = TRUE,
    draws      = control$laplace_draws,
    refresh    = 0L
  )
  if (!is.null(seed)) laplace_args$seed <- as.integer(seed) + 1000L
  lap <- tryCatch(
    do.call(model$laplace, laplace_args),
    error = function(e) NULL
  )
  if (is.null(lap)) {
    gdpar_abort(
      paste(
        "Path C: cmdstanr::laplace() failed on the best init returned by",
        "multi_start optimize(); Hessian likely singular or",
        "non-positive-definite at the candidate joint MAP."
      ),
      class = "gdpar_eb_numerical_error",
      data = list(history_lp = history_lp, best_idx = best_idx)
    )
  }
  draws <- lap$draws()
  all_vars <- dimnames(draws)$variable %||% character(0L)
  # Enumerate the J*K*p theta_ref_kp[g, k, c] variables expected.
  vars_table <- expand.grid(g = seq_len(J), k = seq_len(K),
                             c = seq_len(p), KEEP.OUT.ATTRS = FALSE)
  vars_kp <- sprintf("theta_ref_kp[%d,%d,%d]",
                      vars_table$g, vars_table$k, vars_table$c)
  available <- intersect(vars_kp, all_vars)
  if (length(available) != length(vars_kp)) {
    gdpar_abort(
      sprintf(
        "Path C EB Laplace draws missing theta_ref_kp entries: expected %d found %d.",
        length(vars_kp), length(available)
      ),
      class = "gdpar_internal_error"
    )
  }
  draws_kp <- posterior::subset_draws(draws, variable = vars_kp)
  mat_kp <- posterior::as_draws_matrix(draws_kp)
  hat_kp <- colMeans(mat_kp)
  # Build 3D array [J, K, p] in (g, k, c) row-major order matching the
  # vars_table enumeration.
  theta_ref_kp_hat <- array(NA_real_, dim = c(J, K, p))
  for (i in seq_along(vars_kp)) {
    theta_ref_kp_hat[vars_table$g[i], vars_table$k[i],
                      vars_table$c[i]] <- hat_kp[i]
  }
  # D43 = (a): joint Laplace + block-diagonal extraction per slot. For
  # the per-slot p by p covariance we average over groups when J > 1
  # (canonical aggregation matching the per-slot interpretation of
  # Proposition 7B* of v07b Section 5.1). The structure preserves
  # cross-coord terms within each slot (D37 = (i)) and drops cross-slot
  # terms (block-diagonal canonical assumption under factorized prior).
  laplace_result_per_slot <- vector("list", K)
  names(laplace_result_per_slot) <- paste0("slot_", seq_len(K))
  kappa_per_slot     <- rep(NA_real_, K)
  lm_lambda_per_slot <- rep(0.0, K)
  lm_n_iter_per_slot <- rep(0L, K)
  lm_status_per_slot <- rep("not_needed", K)
  for (k in seq_len(K)) {
    # Per-slot block: rows of vars_table with this k, all groups and
    # all coords.
    slot_rows <- which(vars_table$k == k)
    cov_slot <- if (length(slot_rows) > 1L) {
      stats::cov(mat_kp[, slot_rows, drop = FALSE])
    } else {
      matrix(stats::var(mat_kp[, slot_rows[1L]]), 1L, 1L)
    }
    if (J > 1L) {
      # Average the J group blocks of size p by p into a single block,
      # matching the per-slot canonical aggregation of v07b Section 5.1.
      acc <- matrix(0.0, p, p)
      for (g in seq_len(J)) {
        idx <- which(vars_table$g == g & vars_table$k == k)
        acc <- acc + cov_slot[seq_along(idx) + (g - 1L) * p,
                                seq_along(idx) + (g - 1L) * p,
                                drop = FALSE]
      }
      cov_slot_p <- acc / J
    } else {
      cov_slot_p <- cov_slot
    }
    # O5*-EBFB anti-fragility (Sub-bloque 9.3.b, B9.2): apply the
    # canonical adaptive Levenberg-Marquardt helper per-slot, matching
    # the per-slot canonical aggregation already established by D43 of
    # Sub-fase 8.6.D. Each slot may converge with its own lambda or
    # exhaust independently; the post-ridge kappa is what gets gated
    # below.
    lm_slot <- .gdpar_eb_lm_perturb(cov_slot_p, control)
    cov_slot_p             <- lm_slot$cov_perturbed
    kappa_per_slot[k]      <- lm_slot$kappa_post
    lm_lambda_per_slot[k]  <- lm_slot$lambda_used
    lm_n_iter_per_slot[k]  <- lm_slot$n_iter
    lm_status_per_slot[k]  <- lm_slot$status
    laplace_result_per_slot[[k]] <- list(theta_ref_cov_k = cov_slot_p)
  }
  bad_threshold <- (is.finite(kappa_per_slot) &
                    kappa_per_slot > control$kappa_threshold)
  bad_exhausted <- (lm_status_per_slot == "exhausted")
  if (any(bad_threshold) || any(bad_exhausted)) {
    bad_k <- which(bad_threshold | bad_exhausted)
    gdpar_abort(
      sprintf(
        paste(
          "Path C: slot(s) %s have marginal Hessian condition number kappa",
          "= %s exceeding threshold %.3g (laplace_control$kappa_threshold)",
          "after adaptive Levenberg-Marquardt ridge (status = [%s], lambda",
          "= [%s], n_iter = [%s]). Multivariate joint Laplace is unreliable",
          "for this K x p configuration; consider gdpar() (FB) instead."
        ),
        paste(bad_k, collapse = ", "),
        paste(format(kappa_per_slot[bad_k], digits = 3), collapse = ", "),
        control$kappa_threshold,
        paste(lm_status_per_slot[bad_k], collapse = ", "),
        paste(format(lm_lambda_per_slot[bad_k], digits = 3),
              collapse = ", "),
        paste(lm_n_iter_per_slot[bad_k], collapse = ", ")
      ),
      class = "gdpar_eb_numerical_error",
      data = list(kappa_per_slot     = kappa_per_slot,
                  lm_lambda_per_slot = lm_lambda_per_slot,
                  lm_n_iter_per_slot = lm_n_iter_per_slot,
                  lm_status_per_slot = lm_status_per_slot,
                  history_lp = history_lp)
    )
  }
  history_clean <- history_lp[is.finite(history_lp)]
  dispersion <- if (length(history_clean) > 1L) {
    sd(history_clean) / max(abs(mean(history_clean)), 1)
  } else {
    NA_real_
  }
  if (isTRUE(dispersion > 0.05) && verbose) {
    gdpar_warn(
      sprintf(
        paste(
          "Path C multi-start dispersion across %d Laplace inits is %.3f",
          "(threshold 0.05): possible multimodality of the joint marginal",
          "likelihood under K x p (v07b Section 9.5 O5*-EBFB).",
          "Consider gdpar() (FB)."
        ),
        sum(is.finite(history_lp)), dispersion
      ),
      class = "gdpar_diagnostic_warning",
      data = list(dispersion = dispersion, history_lp = history_lp)
    )
  }
  list(
    theta_ref_kp_hat        = theta_ref_kp_hat,
    laplace_result_per_slot = laplace_result_per_slot,
    diagnostics = list(
      kappa_per_slot              = kappa_per_slot,
      lm_lambda_per_slot          = lm_lambda_per_slot,
      lm_n_iter_per_slot          = lm_n_iter_per_slot,
      lm_status_per_slot          = lm_status_per_slot,
      multi_start_dispersion      = dispersion,
      marginal_log_lik_history    = history_lp,
      best_init_index             = best_idx
    )
  )
}

#' Internal: resolve the anchor argument for the Path C K x p regime
#'
#' Sub-phase 8.6.D (Session 13c 2026-05-26). Mirrors
#' \code{resolve_anchor_K} and \code{resolve_anchor_multi} for the K x p
#' regime: accepts a numeric scalar (broadcast to a K by p matrix), a
#' numeric matrix of shape K by p, or one of \code{"prior_mean"}
#' (zeros) and \code{"empirical_y"} (slot 1 link applied to colMeans of
#' the outcome; slots 2..K anchored at zero).
#'
#' @return Numeric matrix of shape K by p.
#' @keywords internal
#' @noRd
.resolve_anchor_KxP <- function(anchor, family, y_matrix, K, p, verbose) {
  if (is.numeric(anchor) && length(anchor) == 1L && is.finite(anchor)) {
    return(matrix(as.double(anchor), nrow = K, ncol = p))
  }
  if (is.matrix(anchor) && nrow(anchor) == K && ncol(anchor) == p) {
    if (any(!is.finite(anchor))) {
      gdpar_abort(
        "Argument 'anchor' as a matrix must contain only finite values.",
        class = "gdpar_input_error"
      )
    }
    return(matrix(as.double(anchor), K, p))
  }
  if (is.character(anchor) && length(anchor) == 1L) {
    if (anchor == "prior_mean") {
      return(matrix(0.0, K, p))
    }
    if (anchor == "empirical_y") {
      yb <- colMeans(as.matrix(y_matrix))
      eta_loc <- vapply(seq_len(p), function(c) {
        tryCatch(
          family$linkfun(yb[c]),
          error = function(e) {
            gdpar_abort(
              sprintf(
                paste(
                  "Failed to apply the location-slot link to the column",
                  "%d mean of the n by p outcome: %s"
                ),
                c, conditionMessage(e)
              ),
              class = "gdpar_input_error"
            )
          }
        )
      }, numeric(1L))
      out <- matrix(0.0, K, p)
      out[1L, ] <- eta_loc
      if (verbose) {
        gdpar_inform(
          sprintf(
            paste(
              "Path C: empirical anchor on slot 1 = link(colMeans(y)) =",
              "[%s]; every other slot anchored at 0."
            ),
            paste(format(eta_loc, digits = 4), collapse = ", ")
          ),
          class = "gdpar_anchor_message"
        )
      }
      return(out)
    }
  }
  gdpar_abort(
    paste(
      "Argument 'anchor' must be a numeric scalar (broadcast to K x p),",
      "a numeric matrix of shape K x p, 'prior_mean', or 'empirical_y'."
    ),
    class = "gdpar_input_error",
    data = list(received = anchor, K = K, p = p)
  )
}

#' Internal: Path C EB orchestrator (K >= 2 AND p >= 2)
#'
#' Sub-phase 8.6.D (decisions D36 = alpha + D37 = (i) canonized
#' Session 13a 2026-05-25; decisions D38'' = (h) + D39 + D40' + D41
#' canonized Session 13b 2026-05-25; decision D43 = (a) canonized
#' Session 13c 2026-05-26 implementation). Companion of
#' \code{.gdpar_eb_run_K} for the full multivariate K x p regime.
#' Routes the EB pipeline through the dedicated templates
#' \code{inst/stan/amm_eb_marginal_KxP.stan} and
#' \code{inst/stan/amm_eb_conditional_KxP.stan}.
#'
#' Path C consumes a named list of K amm_spec objects each with p >= 2
#' (homogeneous p across slots; heterogeneous-p deferred to Block 9.x)
#' plus a single outcome matrix-column shared across the K slots (per
#' D38'' = (h) canonical composition of Theorem 7C* + Theorem 7A* of
#' v07b). Coverage initial: stan_id in {1 (Gaussian K=2), 3 (NB K=2)};
#' remaining Path B set deferred per D40' + opening Section 6.1
#' numerical caveat.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_run_KxP <- function(amm_list_canonical, family, data, prior,
                              anchor, outcome_name, formula_env,
                              family_id_k_vector,
                              skip_id_check, chains, iter_warmup,
                              iter_sampling, adapt_delta,
                              max_treedepth, refresh, verbose, seed,
                              group, parametrization,
                              id_check_rigor, eb_correction,
                              laplace_control, call, ...) {
  slot_names <- names(amm_list_canonical)
  K <- length(amm_list_canonical)
  p_per_slot <- vapply(amm_list_canonical,
                       function(a) as.integer(a$p %||% 1L),
                       integer(1L))
  p_resolved <- max(p_per_slot)
  if (any(p_per_slot != p_resolved)) {
    gdpar_abort(
      sprintf(
        paste(
          "Path C (Sub-phase 8.6.D first iteration) requires homogeneous p",
          "across the K = %d slots; received p_per_slot = %s.",
          "Heterogeneous p across slots is registered as a follow-on debt",
          "against Block 9.x."
        ),
        K, paste(p_per_slot, collapse = ", ")
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(K = K, p_per_slot = p_per_slot)
    )
  }
  if (p_resolved < 2L) {
    gdpar_abort(
      sprintf(
        "Internal error: Path C dispatcher reached .gdpar_eb_run_KxP with p = %d.",
        p_resolved
      ),
      class = "gdpar_internal_error"
    )
  }
  .gdpar_eb_check_stan_id_for_path(family, K, p_resolved)
  if (!outcome_name %in% colnames(data)) {
    gdpar_abort(
      sprintf("Outcome variable '%s' not found in data.", outcome_name),
      class = "gdpar_input_error"
    )
  }
  y_obj <- data[[outcome_name]]
  if (!is.matrix(y_obj)) {
    gdpar_abort(
      sprintf(
        paste(
          "Path C (K > 1 AND p > 1) requires outcome '%s' to be a matrix",
          "column in 'data' with %d columns shared across the K slots",
          "(D38'' = (h) composition of Theorem 7C* + 7A* of v07b)."
        ),
        outcome_name, p_resolved
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(y_obj), outcome_name = outcome_name)
    )
  }
  if (ncol(y_obj) != p_resolved) {
    gdpar_abort(
      sprintf(
        paste(
          "Path C outcome matrix-column '%s' has %d columns but the K",
          "slots declare p = %d; the outcome must carry exactly p columns",
          "shared across the K slots (D38'' = (h))."
        ),
        outcome_name, ncol(y_obj), p_resolved
      ),
      class = "gdpar_input_error",
      data = list(received_ncol = ncol(y_obj), p_resolved = p_resolved)
    )
  }
  y_obj_bad <- if (is.numeric(y_obj)) !is.finite(y_obj) else is.na(y_obj)
  if (any(y_obj_bad)) {
    gdpar_abort(
      sprintf(
        "Outcome '%s' contains %d non-finite value(s) (NA, NaN or Inf).",
        outcome_name, sum(y_obj_bad)
      ),
      class = "gdpar_input_error"
    )
  }
  y_matrix <- as.matrix(y_obj)
  union_vars <- character(0L)
  for (k in seq_len(K)) {
    a_k <- amm_list_canonical[[k]]
    if (!is.null(a_k$a)) union_vars <- c(union_vars, all.vars(a_k$a))
    if (!is.null(a_k$b)) union_vars <- c(union_vars, all.vars(a_k$b))
    if (!is.null(a_k$dims)) {
      for (d in a_k$dims) {
        if (!is.null(d$a)) union_vars <- c(union_vars, all.vars(d$a))
        if (!is.null(d$b)) union_vars <- c(union_vars, all.vars(d$b))
      }
    }
  }
  union_vars <- unique(union_vars)
  rhs_str <- if (length(union_vars) == 0L) "1" else
             paste(union_vars, collapse = " + ")
  formula_full <- stats::as.formula(
    paste(outcome_name, "~", rhs_str),
    env = formula_env
  )
  rhs <- formula_full[c(1L, 3L)]
  rhs <- stats::update(rhs, ~ . + 0)
  any_W_input <- any(vapply(amm_list_canonical,
                             function(a) !is.null(a$W),
                             logical(1L)))
  if (isTRUE(any_W_input)) {
    gdpar_abort(
      paste(
        "Path C (Sub-phase 8.6.D first iteration) has the modulating W",
        "component disabled per D39. At least one slot declares",
        "W != NULL; remove the W() wrapper or defer to Block 9.x."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  design_KxP <- .build_amm_design_KxP(amm_list_canonical, data,
                                       formula_rhs = rhs)
  anchor_value <- .resolve_anchor_KxP(anchor, family, y_matrix,
                                       K, p_resolved, verbose)
  # Identifiability checks per slot + per coord. The K-level + per-slot
  # multivariate identifiability machinery is non-trivial under Path C
  # and the canonical contract of Charter Section 3.4 lets the first
  # iteration skip ID checks under documented caveat; users can opt in
  # via skip_id_check = FALSE on subsequent iterations once
  # check_identifiability_KxP exists.
  id_report <- if (skip_id_check) {
    NULL
  } else {
    if (verbose) {
      gdpar_inform(
        paste(
          "Path C first iteration: per-slot per-coord identifiability",
          "check is delegated to per-slot multivariate check; the K-level",
          "joint check is deferred to a later iteration of 8.6.D."
        ),
        class = "gdpar_eb_message"
      )
    }
    per_slot_reports <- vector("list", K)
    names(per_slot_reports) <- slot_names
    for (k in seq_len(K)) {
      a_k <- amm_list_canonical[[k]]
      slot_vars <- unique(union(
        if (!is.null(a_k$a)) all.vars(a_k$a) else character(0L),
        if (!is.null(a_k$b)) all.vars(a_k$b) else character(0L)
      ))
      if (!is.null(a_k$dims)) {
        for (d in a_k$dims) {
          if (!is.null(d$a)) slot_vars <- union(slot_vars, all.vars(d$a))
          if (!is.null(d$b)) slot_vars <- union(slot_vars, all.vars(d$b))
        }
      }
      slot_rhs_str <- if (length(slot_vars) == 0L) "1" else
                       paste(slot_vars, collapse = " + ")
      slot_rhs <- stats::as.formula(
        paste("~", slot_rhs_str), env = formula_env
      )
      diag_theta <- anchor_value[k, 1L]
      if (abs(diag_theta) < 1e-8 && !is.null(a_k$b)) diag_theta <- 1
      rep_k <- tryCatch(
        gdpar_check_identifiability(
          amm = a_k, data = data,
          theta_ref_init = diag_theta,
          formula_rhs = slot_rhs
        ),
        error = function(e) {
          if (verbose) {
            gdpar_inform(
              sprintf(
                "Path C ID check skipped for slot '%s' due to: %s",
                slot_names[k], conditionMessage(e)
              ),
              class = "gdpar_eb_message"
            )
          }
          NULL
        }
      )
      per_slot_reports[[k]] <- rep_k
    }
    per_slot_reports
  }
  group_info <- .resolve_group_argument(group, data, n = nrow(y_matrix),
                                         verbose = verbose)
  stan_data <- .assemble_stan_data_KxP(
    design_KxP, family, amm_list_canonical, y_matrix, anchor_value,
    group_id = if (is.null(group_info)) NULL else group_info$group_id
  )
  cp_a <- identical(parametrization, "cp")
  cp_W <- identical(parametrization, "cp")
  parametrization_resolved <- list(
    cp_a       = isTRUE(cp_a),
    cp_W       = isTRUE(cp_W),
    cp_a_per_K = NULL,
    meta = list(
      mode = "eb_KxP_path_C",
      note = paste(
        "Sub-phase 8.6.D Path C first iteration ships with NCP per slot",
        "per coord hardcoded in amm_eb_*_KxP.stan; per-slot preflight",
        "(cp_a_per_K) is queued for Block 9.x."
      ),
      requested = list(parametrization = parametrization)
    )
  )
  marginal_src <- .gdpar_eb_generate_stan_marginal(
    prior,
    cp_a   = parametrization_resolved$cp_a,
    cp_W   = parametrization_resolved$cp_W,
    K      = K,
    p      = p_resolved,
    family = family
  )
  marginal_path <- write_stan_to_tempfile(marginal_src)
  require_suggested("cmdstanr",
                    "fit Path 1 Empirical-Bayes models via cmdstanr::laplace()")
  marginal_model <- cmdstanr::cmdstan_model(marginal_path)
  laplace_result <- .gdpar_eb_maximize_marginal_KxP(
    model     = marginal_model,
    stan_data = stan_data,
    control   = laplace_control,
    seed      = seed,
    verbose   = verbose
  )
  conditional_src <- .gdpar_eb_generate_stan_conditional(
    prior,
    cp_a   = parametrization_resolved$cp_a,
    cp_W   = parametrization_resolved$cp_W,
    K      = K,
    p      = p_resolved,
    family = family
  )
  conditional_path  <- write_stan_to_tempfile(conditional_src)
  conditional_model <- cmdstanr::cmdstan_model(conditional_path)
  stan_data_cond <- stan_data
  # theta_ref_kp_data: cmdstanr's automatic packing accepts an R 3D
  # array of shape [J_groups, K, p] for the Stan declaration
  # array[J_groups, K] vector[p] theta_ref_kp_data.
  stan_data_cond$theta_ref_kp_data <- laplace_result$theta_ref_kp_hat
  sample_args <- list(
    data            = stan_data_cond,
    chains          = as.integer(chains),
    iter_warmup     = as.integer(iter_warmup),
    iter_sampling   = as.integer(iter_sampling),
    adapt_delta     = adapt_delta,
    max_treedepth   = as.integer(max_treedepth),
    refresh         = as.integer(refresh),
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) sample_args$seed <- as.integer(seed)
  extra <- list(...)
  for (nm in names(extra)) sample_args[[nm]] <- extra[[nm]]
  fit_cond <- do.call(conditional_model$sample, sample_args)
  diagnostics <- compute_diagnostics(fit_cond, verbose = verbose)
  correction <- .gdpar_eb_correction_tensor(
    eb_correction           = eb_correction,
    laplace_result_per_slot = laplace_result$laplace_result_per_slot,
    K                       = K,
    p                       = p_resolved,
    verbose                 = verbose
  )
  # Per-group, per-slot, per-coord SE derived from the diagonal of each
  # slot's covariance block (replicated across groups under the
  # per-slot canonical aggregation D43 = (a)).
  theta_ref_kp_se <- array(NA_real_, dim = c(stan_data$J_groups, K,
                                              p_resolved))
  for (k in seq_len(K)) {
    cov_k <- laplace_result$laplace_result_per_slot[[k]]$theta_ref_cov_k
    se_k <- if (is.matrix(cov_k))
              sqrt(pmax(diag(cov_k), 0)) else
              rep(NA_real_, p_resolved)
    for (g in seq_len(stan_data$J_groups)) {
      theta_ref_kp_se[g, k, ] <- se_k
    }
  }
  obj <- list(
    theta_ref_kp_hat              = laplace_result$theta_ref_kp_hat,
    theta_ref_kp_se               = theta_ref_kp_se,
    theta_ref_kp_cov_per_slot     = lapply(
      laplace_result$laplace_result_per_slot,
      function(x) x$theta_ref_cov_k
    ),
    conditional_fit               = fit_cond,
    amm_list_canonical            = amm_list_canonical,
    family                        = family,
    prior                         = prior,
    design_KxP                    = design_KxP,
    anchor                        = anchor_value,
    stan_data                     = stan_data,
    identifiability_report        = id_report,
    diagnostics                   = diagnostics,
    diagnostics_numerical         = laplace_result$diagnostics,
    parametrization               = parametrization_resolved,
    group_info                    = group_info,
    correction_applied            = correction$applied,
    correction_tensor_constant    = correction$constant,
    correction_tensor_dispositions = correction$slot_dispositions,
    call                          = call,
    path                          = "eb_KxP",
    K                             = K,
    p                             = p_resolved,
    slot_names                    = slot_names
  )
  class(obj) <- c("gdpar_eb_fit", "list")
  obj
}

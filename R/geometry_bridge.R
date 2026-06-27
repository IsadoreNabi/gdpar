# Integration of the geometry-adaptive orchestrator with the package fit
# engine (Block RG, RG.6 part ii). Two opt-in, net-new layers; the default
# fit path (gdpar()) is untouched, so the goldens stay bit-identical.
#
#   * gdpar_geom_bridge(object): the durable, path-agnostic core. It takes an
#     already-fitted gdpar_fit, reads its compiled cmdstan model and Stan data,
#     enables the standalone log_prob / grad_log_prob / hessian methods, derives
#     the unconstrained dimension and a posterior-mean warm-start, and returns
#     the (target, geom_target, fisher, reference) tuple that
#     gdpar_geom_orchestrate() consumes. It NEVER touches gdpar(). This is the
#     tool RG.7 will point at the real Tweedie count.
#
#   * gdpar_geom_fit(formula, ...): the one-call ergonomic entry, a SISTER of
#     gdpar() (not an internal branch). It builds and compiles the K-individual
#     model through the SHARED seam .gdpar_K_build() -- the same single source
#     gdpar()'s .gdpar_K() uses, so there is no duplication and no throwaway
#     computation -- assembles the bridge, runs the orchestrator and returns a
#     rich object. Its scope in RG.6 is the K-individual path (where the Tweedie
#     lives); the multi / main paths can be wired the same way later.
#
# Design synthesis decided by the user (session B9.30, canonised under D96/D97):
# the bridge combines the best of the three options once weighed (a thin wrapper,
# a geometry = "auto" branch inside gdpar(), and a standalone bridge) without
# their weaknesses -- the bridge is the durable, zero-coupling core; the wrapper
# is the ergonomic single call; neither freezes a contract that RG.7 might have
# to redesign, and the default branch stays bit-identical.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Enable the standalone gradient/Hessian methods on a cmdstan fit. Idempotent:
# init_model_methods() is a no-op when the methods are already attached. The
# Hessian is requested so the Riemannian / SoftAbs level has $hessian, but its
# compilation is BEST-EFFORT: higher-order autodiff may fail to compile for
# models with custom densities (the Tweedie lpdf is the case RG.7 targets). On
# Hessian failure we fall back to gradient-only methods, which fully serve the
# Euclidean, dense and sub-Riemannian (simulated-Fisher) levels -- exactly the
# levels the count needs. Only a failure to expose log_prob / grad_log_prob at
# all is fatal. Returns whether the Hessian method is available.
.gdpar_geom_enable_methods <- function(fit, hessian = TRUE, seed = 1L) {
  try_init <- function(hess) tryCatch({
    fit$init_model_methods(seed = as.integer(seed), verbose = FALSE,
                           hessian = hess)
    TRUE
  }, error = function(e) e)
  r <- try_init(isTRUE(hessian))
  if (isTRUE(r)) return(invisible(list(has_hessian = isTRUE(hessian))))
  if (isTRUE(hessian)) {
    r2 <- try_init(FALSE)
    if (isTRUE(r2)) {
      gdpar_warn(paste(
        "The standalone Hessian could not be compiled (higher-order autodiff",
        "is unavailable for this model); falling back to gradient-only",
        "methods. The Riemannian SoftAbs level is disabled, but the Euclidean,",
        "dense and sub-Riemannian (expected-Fisher) levels remain available."),
        class = "gdpar_geometry_warning")
      return(invisible(list(has_hessian = FALSE)))
    }
    r <- r2
  }
  gdpar_abort(
    sprintf(paste(
      "Could not expose the standalone log_prob / grad_log_prob methods on the",
      "cmdstan model (init_model_methods failed: %s). The geometry engine needs",
      "them; a C++ toolchain capable of compiling the model methods is",
      "required."), conditionMessage(r)),
    class = "gdpar_input_error")
}

# The unconstrained dimension and posterior-mean warm-start of a fit with
# methods enabled (the unconstrained scale is where the engine integrates).
.gdpar_geom_unconstrained_summary <- function(fit) {
  um <- tryCatch(
    posterior::as_draws_matrix(fit$unconstrain_draws()),
    error = function(e) gdpar_abort(
      sprintf("Failed to unconstrain the posterior draws: %s.",
              conditionMessage(e)),
      class = "gdpar_input_error"))
  um <- as.matrix(um)
  list(dim = ncol(um), reference = unname(colMeans(um)))
}

# Assemble the bridge object from a methods-enabled engine source.
#   model        : a CmdStanModel compiled with methods (the re-samplable target
#                  the diagnostic needs).
#   stan_data    : the Stan data list.
#   dim          : the unconstrained dimension.
#   engine_fit   : a cmdstan fit (with methods) to back the engine geom_target;
#                  when NULL the geom_target is derived from `model` (a cheap
#                  one-iteration sample exposes the methods).
.gdpar_geom_bridge_core <- function(model, stan_data, dim, fisher = NULL,
                                    reference = NULL, engine_fit = NULL,
                                    has_hessian = TRUE, extra = list()) {
  if (!is.null(fisher) && !is.function(fisher)) {
    gdpar_abort("'fisher' must be a function of theta or NULL.",
                class = "gdpar_input_error")
  }
  target <- list(model = model, dim = as.integer(dim), data = stan_data)
  src <- if (is.null(engine_fit)) model else engine_fit
  geom_target <- gdpar_geom_target(object = src, dim = as.integer(dim),
                                   data = stan_data)
  # When the model Hessian was not compiled, drop the hessian closure so the
  # SoftAbs level degrades gracefully instead of erroring at call time.
  if (!isTRUE(has_hessian)) geom_target$hessian <- NULL
  obj <- list(
    target      = target,
    geom_target = geom_target,
    fisher      = fisher,
    reference   = reference,
    dim         = as.integer(dim),
    model       = model,
    stan_data   = stan_data
  )
  obj <- c(obj, extra)
  class(obj) <- c("gdpar_geom_bridge", "list")
  obj
}

#' Bridge a fitted gdpar model to the geometry-adaptive orchestrator
#'
#' Turns an already-fitted \code{\link{gdpar}} object into the inputs the
#' geometry-adaptive controller \code{\link{gdpar_geom_orchestrate}} consumes,
#' \emph{without touching the fit path}. This is the durable, path-agnostic core
#' of the Block RG integration (RG.6 part ii): it reads the compiled cmdstan
#' model and the Stan data carried by the fit, exposes the standalone
#' \code{log_prob} / \code{grad_log_prob} / \code{hessian} methods, derives the
#' unconstrained dimension and a posterior-mean warm-start, and packages a
#' \code{target} (a re-samplable cmdstan model for the size-invariant
#' diagnostic) together with a \code{geom_target} (the engine sampling target on
#' the unconstrained scale).
#'
#' @details
#' The orchestrator needs two things: a target it can \emph{re-sample} for the
#' diagnostic pilots, and an engine target exposing the unconstrained
#' log-density and its gradient (and Hessian, for the Riemannian / SoftAbs
#' level). A fitted \code{gdpar} object carries the posterior draws but its
#' \code{CmdStanMCMC} object cannot be re-sampled. The bridge therefore
#' recompiles a fresh \code{CmdStanModel} from the fit's own Stan source
#' (\code{$code()}) with \code{compile_model_methods = TRUE} (cmdstanr's
#' content-hash cache makes this a cache hit when the methods variant already
#' exists), which serves both consumers, and reads the unconstrained dimension
#' and posterior mean from the fit's draws.
#'
#' The bridge is \strong{path-agnostic}: it works for any \code{gdpar_fit} that
#' carries \code{$fit} (a \code{CmdStanMCMC}) and \code{$stan_data} -- the
#' K-individual, multi-coordinate and single-coordinate paths alike. It is the
#' tool RG.7 points at the real Tweedie count of benchmark 9.2.O.
#'
#' Nothing here modifies \code{\link{gdpar}}; the returned object is plain data
#' plus closures, so the default fit branch stays bit-identical and the goldens
#' are untouched.
#'
#' @param object A fitted \code{gdpar_fit} (the result of \code{\link{gdpar}}).
#' @param fisher Optional function of the unconstrained \code{theta} returning
#'   the expected Fisher information; required by the sub-Riemannian level (the
#'   Tweedie remedy). Use \code{\link{gdpar_geom_fisher_simulator}} when there is
#'   no closed form.
#' @param reference Optional unconstrained reference position (warm-start for the
#'   position-dependent levels). Defaults to the posterior mean read from the
#'   fit.
#' @param hessian Logical; whether to compile the standalone Hessian method
#'   (needed by the Riemannian SoftAbs level). Defaults to \code{TRUE}.
#' @param methods_seed Integer seed forwarded to
#'   \code{init_model_methods()} (the standalone methods are deterministic; this
#'   only seeds any internal RNG). Defaults to \code{1L}.
#' @param ... Reserved for future extension; currently unused.
#'
#' @return An object of class \code{gdpar_geom_bridge}: a list with
#'   \code{target} (the cmdstan-model diagnostic target), \code{geom_target}
#'   (the engine target), \code{fisher}, \code{reference}, \code{dim} (the
#'   unconstrained dimension), \code{model} (the methods-enabled
#'   \code{CmdStanModel}) and \code{stan_data}. Feed these to
#'   \code{\link{gdpar_geom_orchestrate}}.
#'
#' @seealso \code{\link{gdpar_geom_orchestrate}}, \code{\link{gdpar_geom_fit}},
#'   \code{\link{gdpar_geom_target}}, \code{\link{gdpar_geom_fisher_simulator}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 80
#'   x <- rnorm(n); z <- rnorm(n)
#'   y <- rnorm(n, 0.5 + 0.8 * (x - mean(x)), exp(-0.2 + 0.4 * (z - mean(z))))
#'   d <- data.frame(y = y, x = x, z = z)
#'   fit <- gdpar(gdpar_bf(y ~ a(x), sigma ~ a(z)), data = d,
#'                family = gdpar_family("gaussian"), chains = 1,
#'                iter_warmup = 200, iter_sampling = 200, refresh = 0,
#'                seed = 1, skip_id_check = TRUE, verbose = FALSE)
#'   bridge <- gdpar_geom_bridge(fit)
#'   bridge$dim
#'   # res <- gdpar_geom_orchestrate(bridge$target, bridge$geom_target,
#'   #                               reference = bridge$reference)
#' }
#' }
#'
#' @export
gdpar_geom_bridge <- function(object, fisher = NULL, reference = NULL,
                              hessian = TRUE, methods_seed = 1L, ...) {
  if (!inherits(object, "gdpar_fit")) {
    gdpar_abort("'object' must be a fitted gdpar model (class 'gdpar_fit').",
                class = "gdpar_input_error")
  }
  csfit <- object$fit
  if (!inherits(csfit, c("CmdStanMCMC", "CmdStanFit"))) {
    gdpar_abort(paste(
      "The fit does not carry a cmdstan MCMC object in $fit; the geometry",
      "bridge needs a Path 1 (Bayesian) gdpar fit."),
      class = "gdpar_input_error")
  }
  stan_data <- object$stan_data
  if (is.null(stan_data)) {
    gdpar_abort("The fit does not carry $stan_data; cannot build the bridge.",
                class = "gdpar_input_error")
  }
  require_suggested("cmdstanr", "bridge a gdpar fit to the geometry engine")
  require_suggested("posterior", "unconstrain the posterior draws")
  if (!is.null(reference) && !is.numeric(reference)) {
    gdpar_abort("'reference' must be a numeric vector or NULL.",
                class = "gdpar_input_error")
  }

  # Enable the standalone methods on the fitted object (reuses its compiled
  # executable) and read the unconstrained dimension + posterior-mean warm-start.
  meth <- .gdpar_geom_enable_methods(csfit, hessian = hessian,
                                     seed = methods_seed)
  us <- .gdpar_geom_unconstrained_summary(csfit)
  d <- us$dim
  if (is.null(reference)) {
    reference <- us$reference
  } else if (length(reference) != d) {
    gdpar_abort(sprintf("'reference' must have length %d (the unconstrained dimension).", d),
                class = "gdpar_input_error")
  }

  # Recompile a re-samplable CmdStanModel (with methods) from the fit's own
  # Stan source for the diagnostic; cmdstanr's hash cache makes this a hit.
  model <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(csfit$code()),
    compile_model_methods = TRUE)

  # The engine target is backed by the fitted object itself (methods already
  # enabled -> no extra sample); the diagnostic target by the re-samplable model.
  .gdpar_geom_bridge_core(model = model, stan_data = stan_data, dim = d,
                          fisher = fisher, reference = reference,
                          engine_fit = csfit, has_hessian = meth$has_hessian)
}

#' Print method for gdpar_geom_bridge objects
#'
#' @param x A \code{gdpar_geom_bridge}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_bridge <- function(x, ...) {
  cat("<gdpar_geom_bridge>\n")
  cat("  unconstrained dim: ", x$dim, "\n", sep = "")
  cat("  fisher supplied:   ", !is.null(x$fisher), "\n", sep = "")
  cat("  reference:         ",
      if (is.null(x$reference)) "zeros (default)" else
        paste0("length ", length(x$reference), " (warm-start)"), "\n", sep = "")
  cat("  feed (target, geom_target, fisher, reference) to gdpar_geom_orchestrate()\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# gdpar_geom_fit(): one-call ergonomic entry (K-individual scope)
# ---------------------------------------------------------------------------

# Resolve the K-individual front-end inputs the way gdpar() does, reusing the
# same helpers, without touching gdpar(). Returns the arguments .gdpar_K_build()
# needs. Scoped to K > 1 (the distributional / Tweedie regime).
.gdpar_geom_fit_resolve_K <- function(formula, family, amm, W, data) {
  family_is_named_list <- is.list(family) &&
    !inherits(family, "gdpar_family") &&
    !inherits(family, "gdpar_family_multi") &&
    !is.null(names(family)) && all(nzchar(names(family))) &&
    anyDuplicated(names(family)) == 0L &&
    all(vapply(family, inherits, logical(1L), "gdpar_family"))
  if (!family_is_named_list) {
    assert_inherits(family, c("gdpar_family", "gdpar_family_multi"), "family")
  }
  formula_set_input <- inherits(formula, "gdpar_formula_set")
  amm_list_input <- is.list(amm) && !inherits(amm, "amm_spec") &&
    !is.null(names(amm))
  classic_with_amm <- !formula_set_input && inherits(formula, "formula") &&
    length(formula) == 3L && .gdpar_rhs_has_amm_calls(formula)
  if (!(formula_set_input || amm_list_input || classic_with_amm)) {
    gdpar_abort(paste(
      "gdpar_geom_fit() targets the K-individual path: supply a gdpar_bf()",
      "formula set, a named list of amm_spec, or a classic formula whose RHS",
      "carries a()/b()/W() wrapper calls. The single-coordinate and",
      "multi-coordinate paths are not yet wrapped (RG.6 scope = K path)."),
      class = "gdpar_input_error")
  }
  if (formula_set_input) {
    if (!.gdpar_is_default_amm_spec(amm)) {
      gdpar_abort(paste(
        "When 'formula' is a gdpar_formula_set, 'amm' must stay at its default",
        "amm_spec()."), class = "gdpar_input_error")
    }
    amm_list_canonical <- .gdpar_formula_set_to_amm_spec_list(formula, W)
    outcome_name <- formula$outcome
    formula_env <- formula$env
  } else if (amm_list_input) {
    for (slot in names(amm)) {
      if (!nzchar(slot) || !inherits(amm[[slot]], "amm_spec")) {
        gdpar_abort(paste(
          "Every entry of the named list passed to 'amm' must be a named",
          "amm_spec."), class = "gdpar_input_error")
      }
    }
    if (!inherits(formula, "formula") || length(formula) != 3L) {
      gdpar_abort(paste(
        "When 'amm' is a named list of amm_spec, 'formula' must be a",
        "two-sided formula 'y ~ ...'."), class = "gdpar_input_error")
    }
    amm_list_canonical <- amm
    outcome_name <- as.character(formula[[2L]])
    formula_env <- environment(formula)
  } else {
    if (!.gdpar_is_default_amm_spec(amm)) {
      gdpar_abort(paste(
        "The RHS of 'formula' contains AMM wrapper calls; 'amm' must stay at",
        "its default amm_spec()."), class = "gdpar_input_error")
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
  if (K <= 1L) {
    gdpar_abort(paste(
      "gdpar_geom_fit() requires K > 1 (a distributional / K-individual",
      "specification with at least two slots). For K = 1 fit with gdpar() and",
      "bridge the result with gdpar_geom_bridge()."),
      class = "gdpar_input_error", data = list(K = K))
  }
  if (family_is_named_list) {
    het <- .gdpar_resolve_heterogeneous_family_K(family, names(amm_list_canonical))
    family_promoted <- het$location_family
    family_id_k_vector <- het$family_id_k_vector
  } else {
    family_promoted <- .gdpar_promote_scope_per_observation(
      family, names(amm_list_canonical))
    family_id_k_vector <- NULL
  }
  list(amm_list_canonical = amm_list_canonical, family = family_promoted,
       outcome_name = outcome_name, formula_env = formula_env,
       family_id_k_vector = family_id_k_vector, K = K)
}

#' One-call geometry-adaptive fit (K-individual path)
#'
#' The ergonomic single-call entry of the Block RG integration: a \emph{sister}
#' of \code{\link{gdpar}} (not an internal branch of it) that builds and compiles
#' a K-individual model, then runs the opt-in geometry-adaptive controller
#' \code{\link{gdpar_geom_orchestrate}} on it instead of the default NUTS fit. It
#' diagnoses the posterior geometry, selects a level of the sampler hierarchy
#' (Euclidean diagonal / dense, Riemannian, relativistic, sub-Riemannian),
#' samples, re-diagnoses and either resolves or emits a certified limit.
#'
#' @details
#' \code{gdpar_geom_fit()} shares the model-building seam \code{.gdpar_K_build()}
#' with \code{\link{gdpar}}'s internal K path: the model and data are assembled
#' and compiled exactly once, through a single source, with no duplication and
#' no throwaway computation. The model is compiled with the standalone gradient
#' and Hessian methods exposed (\code{compile_model_methods = TRUE}) so the
#' geometry engine can integrate on the unconstrained scale; the default
#' \code{\link{gdpar}} branch, which compiles without those methods, is byte for
#' byte unchanged and its goldens stay bit-identical.
#'
#' This entry is scoped to the K-individual (distributional) regime where the
#' Tweedie count lives. For an already-fitted model, or for the single- and
#' multi-coordinate paths, use \code{\link{gdpar_geom_bridge}} on the fit.
#'
#' Correctness versus efficiency (the honesty convention of ORPHEUS-PIMC
#' section 16.3 the package follows): every sampler level is Metropolis-exact, so
#' \emph{which} geometry is selected only governs efficiency, never the validity
#' of the returned draws.
#'
#' @param formula A \code{gdpar_bf()} formula set, a classic formula with AMM
#'   wrapper calls on the RHS, or (with \code{amm}) a two-sided \code{y ~ ...}.
#' @param family A \code{gdpar_family} (broadcast across the K slots) or a named
#'   list of \code{gdpar_family} (heterogeneous families per slot).
#' @param amm A named list of \code{amm_spec} (one per slot) for the named-list
#'   path; otherwise left at its default.
#' @param W Optional modulating basis for the K-individual paths.
#' @param data A data frame.
#' @param prior A \code{gdpar_prior}; defaults to \code{gdpar_prior()}.
#' @param anchor The slot anchor(s); see \code{\link{gdpar}}.
#' @param skip_id_check Logical; skip the identifiability check.
#' @param parametrization,parametrization_a,parametrization_W The CP/NCP toggles
#'   forwarded to the model build.
#' @param id_check_rigor One of \code{"full"} / \code{"fast"}.
#' @param group Optional grouping variable for grouped anchors.
#' @param fisher Optional expected-Fisher function (sub-Riemannian level). Use
#'   \code{\link{gdpar_geom_fisher_simulator}} when there is no closed form.
#' @param budget,criteria The orchestrator budget / success gate; see
#'   \code{\link{gdpar_geom_orchestrate_budget}} and
#'   \code{\link{gdpar_geom_orchestrate_criteria}}.
#' @param entry_level,level_map Optional overrides of the entry level and the
#'   pathology-to-level map.
#' @param reference Optional unconstrained warm-start position.
#' @param speed,rest_mass The relativistic level's speed and rest mass.
#' @param laplace_fallback Logical; forwarded to
#'   \code{\link{gdpar_geom_orchestrate}}. When \code{TRUE} and the run ends in a
#'   certified limit, a \code{\link{gdpar_geom_laplace}} approximation is attached
#'   (\code{$laplace}) and the status becomes \code{"certified_limit_laplace"}.
#'   Defaults to \code{FALSE} (bit-identical output).
#' @param laplace_draws Number of iid Laplace draws carried on the fallback when
#'   \code{laplace_fallback = TRUE} (default \code{0L}).
#' @param n_grid Optional diagnostic size grid (forwarded).
#' @param seed Integer base seed for the (deterministic) adaptive trajectory.
#' @param verbose Logical; opt-in progress trace.
#' @param ... Forwarded to \code{.gdpar_K_build()} / the diagnostic.
#'
#' @return An object of class \code{gdpar_geom_fit}: a list carrying the
#'   \code{orchestration} (a \code{\link{gdpar_geom_orchestrate}} result), the
#'   \code{bridge}, the \code{status}, and, when resolved, the winning
#'   \code{level}, \code{metric} and \code{draws}; when the budget is exhausted,
#'   the \code{certificate} (and, under \code{laplace_fallback = TRUE}, the
#'   \code{laplace} approximation with status \code{"certified_limit_laplace"}).
#'   It also carries \code{stan_data}, \code{family}, \code{K}, \code{slot_names}
#'   and the \code{call}.
#'
#' @seealso \code{\link{gdpar_geom_bridge}}, \code{\link{gdpar_geom_orchestrate}},
#'   \code{\link{gdpar}}.
#'
#' @examples
#' b <- gdpar_geom_orchestrate_budget()
#' b$max_rounds
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 80
#'   x <- rnorm(n); z <- rnorm(n)
#'   y <- rnorm(n, 0.5 + 0.8 * (x - mean(x)), exp(-0.2 + 0.4 * (z - mean(z))))
#'   d <- data.frame(y = y, x = x, z = z)
#'   b$tune_epsilon <- FALSE
#'   b$probe_iter <- 60L; b$full_iter <- 80L; b$full_warmup <- 80L
#'   res <- gdpar_geom_fit(gdpar_bf(y ~ a(x), sigma ~ a(z)), data = d,
#'                         family = gdpar_family("gaussian"),
#'                         skip_id_check = TRUE, budget = b, n_grid = 1,
#'                         verbose = FALSE)
#'   res$status
#' }
#' }
#'
#' @export
gdpar_geom_fit <- function(formula, family = gdpar_family("gaussian"),
                           amm = amm_spec(), W = NULL, data, prior = NULL,
                           anchor = "prior_mean", skip_id_check = FALSE,
                           parametrization = c("auto", "ncp", "cp"),
                           parametrization_a = NULL, parametrization_W = NULL,
                           id_check_rigor = c("full", "fast"), group = NULL,
                           fisher = NULL, budget = NULL, criteria = NULL,
                           entry_level = NULL, level_map = NULL,
                           reference = NULL, speed = 10, rest_mass = 1,
                           laplace_fallback = FALSE, laplace_draws = 0L,
                           n_grid = NULL, seed = 20260603L, verbose = TRUE,
                           ...) {
  call <- match.call()
  parametrization <- match.arg(parametrization)
  id_check_rigor <- match.arg(id_check_rigor)
  if (!is.null(parametrization_a)) {
    parametrization_a <- match.arg(parametrization_a, c("ncp", "cp"))
  }
  if (!is.null(parametrization_W)) {
    parametrization_W <- match.arg(parametrization_W, c("ncp", "cp"))
  }
  assert_data_frame(data, "data")
  require_suggested("cmdstanr", "fit Path 1 hierarchical Bayesian models")
  require_suggested("posterior", "unconstrain the posterior draws")

  res <- .gdpar_geom_fit_resolve_K(formula, family, amm, W, data)

  inform <- function(msg) if (isTRUE(verbose))
    gdpar_inform(msg, class = "gdpar_optin_message")
  inform("gdpar_geom_fit builds the model once, then runs the opt-in geometry-adaptive controller (it never alters gdpar()'s default fit path).")

  bd <- .gdpar_K_build(
    amm_list_canonical = res$amm_list_canonical, family = res$family,
    data = data, prior = prior, anchor = anchor,
    outcome_name = res$outcome_name, formula_env = res$formula_env,
    skip_id_check = skip_id_check, verbose = verbose, group = group,
    parametrization = parametrization, parametrization_a = parametrization_a,
    parametrization_W = parametrization_W, id_check_rigor = id_check_rigor,
    family_id_k_vector = res$family_id_k_vector,
    compile_model_methods = TRUE)

  # One throwaway-free single-iteration sample exposes the methods and yields
  # the unconstrained dimension; the same fit backs the engine geom_target.
  mfit <- bd$cs_model$sample(
    data = bd$stan_data, chains = 1L, iter_warmup = 1L, iter_sampling = 1L,
    refresh = 0L, show_messages = FALSE, show_exceptions = FALSE)
  meth <- .gdpar_geom_enable_methods(mfit, hessian = TRUE)
  d <- .gdpar_geom_unconstrained_summary(mfit)$dim

  bridge <- .gdpar_geom_bridge_core(
    model = bd$cs_model, stan_data = bd$stan_data, dim = d, fisher = fisher,
    reference = reference, engine_fit = mfit, has_hessian = meth$has_hessian)

  orch <- gdpar_geom_orchestrate(
    target = bridge$target, geom_target = bridge$geom_target,
    fisher = bridge$fisher, reference = bridge$reference,
    level_map = level_map, entry_level = entry_level, budget = budget,
    criteria = criteria, speed = speed, rest_mass = rest_mass,
    laplace_fallback = laplace_fallback, laplace_draws = laplace_draws,
    n_grid = n_grid, seed = seed, verbose = verbose)

  obj <- list(
    orchestration = orch,
    bridge        = bridge,
    status        = orch$status,
    level         = orch$level,
    metric        = orch$metric,
    draws         = orch$draws,
    laplace       = orch$laplace,
    certificate   = orch$certificate,
    stan_data     = bd$stan_data,
    family        = res$family,
    K             = res$K,
    slot_names    = bd$meta$slot_names,
    call          = call
  )
  class(obj) <- c("gdpar_geom_fit", "list")
  obj
}

#' Print method for gdpar_geom_fit objects
#'
#' @param x A \code{gdpar_geom_fit}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_fit <- function(x, ...) {
  cat("<gdpar_geom_fit> K = ", x$K, " | status: ", x$status, "\n", sep = "")
  if (identical(x$status, "resolved")) {
    cat("  resolved at level: ", x$level, "\n", sep = "")
    if (!is.null(x$draws)) {
      cat("  draws: ", nrow(x$draws), " x ", ncol(x$draws), "\n", sep = "")
    }
  } else if (!is.null(x$certificate)) {
    cat("  verdict: ", x$certificate$verdict, "\n", sep = "")
    cat("  prescription items: ", length(x$certificate$prescription), "\n",
        sep = "")
  }
  if (!is.null(x$laplace)) {
    cat("  laplace fallback: '", x$laplace$fit_quality_label, "'\n", sep = "")
  }
  invisible(x)
}

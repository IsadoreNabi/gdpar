# Posterior-geometry forensic diagnostic (Block RG, RG.1).
#
# Runs cheap NUTS pilots through cmdstan and reads size-invariant sampler
# signals (divergences, energy / E-BFMI, tree depth, condition number,
# step/scale ratio) -- NEVER R-hat or ESS on short runs, which were the false
# positives of sessions B9.20/B9.21. A difficulty-vs-n curve separates
# quasi-determinism (difficulty grows with n) from structural pathology
# (difficulty constant in n; the 28-GB nuance). A culprit localiser automates
# the manual probe of B9.21. A transparent rule-based classifier maps the
# signals to a pathology class and the geometry level that remedies it. The
# default fit path of the package is untouched: this is an opt-in, standalone
# forensic tool.

#' Default classifier thresholds for the posterior-geometry diagnostic
#'
#' The decision thresholds of \code{\link{gdpar_geometry_diagnostic}}'s
#' rule-based classifier. They are exposed as data (not hard-coded) so the
#' calibration of RG.1 can tune them against \code{\link{gdpar_geometry_suite}}
#' and so that a user can re-calibrate them for their own setting.
#'
#' @details
#' These defaults were calibrated in sub-phase RG.1.c (session B9.23) against
#' the synthetic suite over a difficulty x pilot-budget x replica grid, with the
#' thresholds proposed by a data-driven Youden cut, regularised to interpretable
#' values on a calibration fold and validated out-of-sample on a held-out fold
#' (held-out balanced accuracy rose from 0.60 to 0.89). The chief change versus
#' the initial hand-set values reflects that the SHORT pilots used by the
#' diagnostic attenuate the signals relative to their asymptotic values (the
#' empirical condition number underestimates the true one, sample kurtosis is
#' damped, boundary pile-up is subtle), so several cuts moved to catch the
#' attenuated signal: \code{condition_high} 50 -> 12, \code{heavy_kurtosis_high}
#' 3 -> 1.8, \code{boundary_prox_high} 0.10 -> 0.02, \code{nslope_grows}
#' 0.30 -> 0.80, \code{funnel_ebfmi_low} 0.25 -> 0.35, \code{heavy_cond_max}
#' 8 -> 25, \code{flat_var_high} 1000 -> 600, \code{multimodal_high} 2 -> 2.5.
#' The calibration is against an idealised synthetic suite and is not a claim of
#' optimality on real posteriors; the funnel and heavy-tail classes remain
#' mutually confusable (per-class recall ~0.6-0.7), which is reported honestly
#' in \code{inst/benchmarks/results/block_rg_calibration.md}. Re-calibrate for a
#' specific application if needed.
#'
#' @return A named list of numeric thresholds.
#' @seealso \code{\link{gdpar_geometry_diagnostic}}.
#' @examples
#' str(gdpar_geometry_thresholds())
#' @export
gdpar_geometry_thresholds <- function() {
  list(
    divergent_rate_high   = 0.01,
    funnel_ebfmi_low      = 0.35,
    heavy_cond_max        = 25,
    treedepth_sat_high    = 0.20,
    condition_high        = 12,
    step_scale_ratio_low  = 0.10,
    nslope_grows          = 0.80,
    flat_var_high         = 600,
    boundary_prox_high    = 0.02,
    boundary_eps          = 0.01,
    multimodal_high       = 2.5,
    heavy_kurtosis_high   = 1.8,
    target_ess            = 400
  )
}

#' Forensic diagnostic of posterior geometry (opt-in)
#'
#' Probe the geometry of a posterior with cheap Hamiltonian pilots and classify
#' the sampling pathology, localise the culprit parameter(s), estimate the
#' difficulty-vs-n behaviour, and recommend the geometry level that remedies it.
#' Built for the Block RG capability: when the no-U-turn sampler stalls (the
#' eBird count / Tweedie case), this tells you \emph{why} and \emph{what} to
#' escalate to, rather than leaving an unexplained \code{rhat = Inf}.
#'
#' The diagnostic uses only size-invariant sampler signals -- the divergence
#' rate, the minimum energy Bayesian fraction of missing information
#' (\dQuote{E-BFMI}), the tree-depth saturation rate, the posterior condition
#' number, and the adapted-step-to-scale ratio. It deliberately does \emph{not}
#' use R-hat or the effective sample size as decision signals: on the short
#' pilots used here those are unreliable and were the false positives of
#' sessions B9.20/B9.21.
#'
#' @param target The posterior to probe. One of three forms (the three-way
#'   adapter): (i) a \code{gdpar_geometry_target} from
#'   \code{\link{gdpar_geometry_suite}} (carries Stan code, a size knob and a
#'   ground-truth label); (ii) a list \code{list(stan_code, stan_data)} or a
#'   compiled \pkg{cmdstanr} model wrapped as \code{list(model = , data = ,
#'   data_n_fn = )}, where the optional \code{data_n_fn(n)} returns the data
#'   list at size \code{n}; (iii) a list \code{list(type = "gdpar", formula = ,
#'   amm = , data = , ...)} probing a real gdpar specification (data is
#'   subsampled to realise the size knob).
#' @param n_grid Optional numeric vector of size-knob values at which to run
#'   pilots. Defaults to the target's own \code{n_grid} (suite targets) or to a
#'   single pilot (other forms without a size knob). At least two values are
#'   required to estimate the difficulty-vs-n curve.
#' @param difficulty Optional pathology-intensity knob forwarded to a suite
#'   target's \code{make()}. Defaults to the target's
#'   \code{default_difficulty}.
#' @param pilot_warmup,pilot_sampling Integer warmup / sampling iterations per
#'   pilot. Kept small on purpose (defaults 150 / 150).
#' @param chains Integer number of chains per pilot. Defaults to 4 (multiple
#'   chains are needed for the multimodality signal).
#' @param adapt_delta Numeric target acceptance for the pilots. Defaults to
#'   0.8 (the diagnostic measures the geometry as the default sampler sees it).
#' @param max_treedepth Integer maximum tree depth for the pilots. Defaults to
#'   10.
#' @param seed Integer base seed. Pilot \code{i} uses \code{seed + i} for
#'   reproducibility.
#' @param thresholds Named list of classifier thresholds; see
#'   \code{\link{gdpar_geometry_thresholds}}. Defaults to that function's
#'   output.
#' @param verbose Logical; print an opt-in cost message before sampling.
#'   Defaults to TRUE.
#' @param ... Additional arguments forwarded to the underlying sampler / fit.
#'
#' @return A list of class \code{gdpar_geometry_diagnostic} with components
#'   \code{pathology} (classified class), \code{confidence},
#'   \code{recommended_geometry} (the remedy level), \code{signals} (data frame
#'   of per-n size-invariant signals), \code{difficulty_curve} (list with
#'   \code{slope}, \code{grows_with_n}), \code{culprit} (ranked data frame),
#'   \code{cost} (list with the cost extrapolation and tractability verdict),
#'   \code{reproducibility} (seed and pilot configuration), and, when the
#'   target carries a ground truth, \code{ground_truth} and \code{correct}.
#'   A \code{print} method summarises the verdict.
#'
#' @section Pathology taxonomy and remedies:
#' The classifier returns one of: \code{isotropic} (Euclidean diagonal;
#' default is fine), \code{anisotropic} (dense Euclidean metric),
#' \code{funnel} (Riemannian metric), \code{heavy_tails} (Finsler /
#' relativistic kinetic energy), \code{quasi_deterministic} (sub-Riemannian),
#' \code{multimodal} (tempering), \code{boundary} (boundary reparametrisation),
#' or \code{flat_direction} (reparametrise / eliminate; the Option A case). The
#' rules are documented in the package source and are calibrated against
#' \code{\link{gdpar_geometry_suite}}.
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} for the pilots and \pkg{posterior} to read draws and
#' sampler diagnostics.
#'
#' @seealso \code{\link{gdpar_geometry_suite}},
#'   \code{\link{gdpar_geometry_thresholds}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   suite <- gdpar_geometry_suite()
#'   diag <- gdpar_geometry_diagnostic(suite$G2_funnel,
#'                                     pilot_warmup = 150, pilot_sampling = 150)
#'   print(diag)
#' }
#' }
#'
#' @export
gdpar_geometry_diagnostic <- function(target, n_grid = NULL, difficulty = NULL,
                                      pilot_warmup = 150L, pilot_sampling = 150L,
                                      chains = 4L, adapt_delta = 0.8,
                                      max_treedepth = 10L, seed = 20260602L,
                                      thresholds = NULL, verbose = TRUE, ...) {
  assert_count(pilot_warmup, "pilot_warmup")
  assert_count(pilot_sampling, "pilot_sampling")
  assert_count(chains, "chains")
  assert_count(max_treedepth, "max_treedepth")
  assert_numeric_scalar(adapt_delta, "adapt_delta", lower = 0, upper = 1)
  assert_count(seed, "seed")
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort("Argument 'verbose' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  require_suggested("cmdstanr", "run posterior-geometry pilots")
  require_suggested("posterior", "read draws and sampler diagnostics")
  if (is.null(thresholds)) thresholds <- gdpar_geometry_thresholds()

  norm <- .gdpar_geom_normalize_target(target, difficulty)
  if (is.null(n_grid)) n_grid <- norm$n_grid_default
  n_grid <- sort(unique(as.numeric(n_grid)))

  controls <- list(
    warmup = pilot_warmup, sampling = pilot_sampling, chains = chains,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth
  )

  if (verbose) {
    gdpar_inform(
      sprintf(
        paste(
          "gdpar_geometry_diagnostic is an opt-in forensic probe.",
          "It will run %d pilot fit(s) (%d size(s) x %d chains x %d+%d iters).",
          "It does not modify any fit; it classifies the posterior geometry",
          "and recommends a sampler level."
        ),
        length(n_grid), length(n_grid), chains, pilot_warmup, pilot_sampling
      ),
      class = "gdpar_optin_message"
    )
  }

  pilots <- vector("list", length(n_grid))
  signal_rows <- list()
  for (i in seq_along(n_grid)) {
    pr <- .gdpar_geom_run_pilot(norm, n_grid[i], controls, seed + i, ...)
    pilots[[i]] <- pr
    signal_rows[[i]] <- data.frame(
      n = n_grid[i],
      divergent_rate = pr$signals$divergent_rate,
      ebfmi_min = pr$signals$ebfmi_min,
      treedepth_sat_rate = pr$signals$treedepth_sat_rate,
      condition_number = pr$signals$condition_number,
      step_scale_ratio = pr$signals$step_scale_ratio,
      lambda_max_cov = pr$signals$lambda_max_cov,
      mean_leapfrog = pr$signals$mean_leapfrog,
      multimodality = pr$signals$multimodality,
      heavy_kurtosis = pr$signals$heavy_kurtosis,
      boundary_proximity = pr$signals$boundary_proximity,
      failed = pr$failed
    )
  }
  signals <- do.call(rbind, signal_rows)

  ncurve <- .gdpar_geom_difficulty_curve(signals, thresholds)
  pilot_max <- pilots[[length(pilots)]]
  culprit <- pilot_max$culprit
  # Multimodality is n-invariant; aggregate it across pilots so a single
  # unlucky seed (all chains landing in one mode) does not mask it.
  sig_classify <- pilot_max$signals
  sig_classify$multimodality <- suppressWarnings(
    max(signals$multimodality, na.rm = TRUE))
  cls <- .gdpar_geom_classify(sig_classify, ncurve, culprit, thresholds)
  cost <- .gdpar_geom_cost(pilot_max, ncurve)

  obj <- list(
    pathology = cls$pathology,
    confidence = cls$confidence,
    recommended_geometry = cls$remedy,
    geometry_level = cls$level,
    signals = signals,
    difficulty_curve = ncurve,
    culprit = culprit,
    cost = cost,
    rule_trace = cls$trace,
    reproducibility = list(
      seed = seed, n_grid = n_grid, controls = controls,
      target_id = norm$meta$id,
      gdpar_version = as.character(utils::packageVersion("gdpar")),
      cmdstan_version = tryCatch(cmdstanr::cmdstan_version(),
                                 error = function(e) NA_character_)
    )
  )
  if (!is.null(norm$meta$ground_truth)) {
    obj$ground_truth <- norm$meta$ground_truth
    obj$correct <- identical(cls$pathology, norm$meta$ground_truth$pathology)
  }
  class(obj) <- c("gdpar_geometry_diagnostic", "list")
  obj
}

# --- target normalisation (three-way adapter) ------------------------------
.gdpar_geom_normalize_target <- function(target, difficulty) {
  # Form (i): a suite target carrying a make() and ground truth.
  if (inherits(target, "gdpar_geometry_target")) {
    diff <- if (is.null(difficulty)) target$default_difficulty else difficulty
    return(list(
      kind = "suite",
      n_grid_default = target$n_grid,
      meta = list(
        id = target$id,
        bounds = target$bounds,
        ground_truth = list(
          pathology = target$pathology,
          geometry_remedy = target$geometry_remedy,
          culprit = target$culprit,
          difficulty_scales_with_n = target$difficulty_scales_with_n
        )
      ),
      make = function(n) target$make(n, diff)
    ))
  }
  # Form (iii): a gdpar specification probed via subsampling.
  if (is.list(target) && identical(target$type, "gdpar")) {
    if (is.null(target$data) || is.null(target$formula)) {
      gdpar_abort(
        "A gdpar target must supply at least 'formula' and 'data'.",
        class = "gdpar_input_error"
      )
    }
    n_full <- nrow(target$data)
    return(list(
      kind = "gdpar",
      n_grid_default = unique(round(exp(seq(
        log(max(30, ceiling(n_full / 4))), log(n_full), length.out = 3)))),
      meta = list(id = "gdpar_spec", bounds = NULL, ground_truth = NULL),
      gdpar_target = target, n_full = n_full
    ))
  }
  # Form (ii): a raw Stan target (code + data) or a compiled cmdstanr model.
  if (is.list(target) && !is.null(target$stan_code)) {
    return(list(
      kind = "stan",
      n_grid_default = 1,
      meta = list(id = "stan_target",
                  bounds = target$bounds %||% NULL, ground_truth = NULL),
      stan_code = target$stan_code,
      data_n_fn = target$data_n_fn %||% function(n) target$stan_data
    ))
  }
  if (is.list(target) && (inherits(target$model, "CmdStanModel") ||
                          !is.null(target$model))) {
    return(list(
      kind = "model",
      n_grid_default = 1,
      meta = list(id = "cmdstan_model",
                  bounds = target$bounds %||% NULL, ground_truth = NULL),
      model = target$model,
      data_n_fn = target$data_n_fn %||% function(n) target$data
    ))
  }
  gdpar_abort(
    paste(
      "Unrecognised 'target'. Supply a gdpar_geometry_target, a",
      "list(stan_code, stan_data), a list(model, data), or a",
      "list(type = 'gdpar', formula, amm, data)."
    ),
    class = "gdpar_input_error"
  )
}

# `%||%` is defined canonically in R/preflight_multi.R.

# --- pilot runner ----------------------------------------------------------
.gdpar_geom_run_pilot <- function(norm, n_knob, controls, seed, ...) {
  fit <- tryCatch(
    .gdpar_geom_sample(norm, n_knob, controls, seed, ...),
    error = function(e) {
      gdpar_warn(
        sprintf("Geometry pilot at n = %s failed: %s.",
                format(n_knob), conditionMessage(e)),
        class = "gdpar_diagnostic_warning"
      )
      NULL
    }
  )
  if (is.null(fit)) {
    return(list(failed = TRUE, signals = .gdpar_geom_na_signals(),
                culprit = .gdpar_geom_empty_culprit()))
  }
  bounds <- norm$meta$bounds
  sig <- .gdpar_geom_extract_signals(fit$cs_fit, fit$pm, fit$chain_id,
                                     bounds, controls$max_treedepth, fit$elapsed)
  culprit <- .gdpar_geom_localize(fit$cs_fit, fit$pm, fit$chain_id, bounds)
  list(failed = FALSE, signals = sig, culprit = culprit,
       elapsed = fit$elapsed)
}

.gdpar_geom_sample <- function(norm, n_knob, controls, seed, ...) {
  if (norm$kind == "gdpar") {
    return(.gdpar_geom_sample_gdpar(norm, n_knob, controls, seed, ...))
  }
  if (norm$kind == "suite") {
    inst <- norm$make(n_knob)
    code <- inst$stan_code
    data <- inst$stan_data
    mod <- .gdpar_geom_compile(code)
  } else if (norm$kind == "stan") {
    code <- norm$stan_code
    data <- norm$data_n_fn(n_knob)
    mod <- .gdpar_geom_compile(code)
  } else { # model
    mod <- norm$model
    data <- norm$data_n_fn(n_knob)
  }
  t0 <- Sys.time()
  cs_fit <- mod$sample(
    data = data, chains = controls$chains,
    parallel_chains = controls$chains,
    iter_warmup = controls$warmup, iter_sampling = controls$sampling,
    adapt_delta = controls$adapt_delta, max_treedepth = controls$max_treedepth,
    seed = seed, refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
    ...
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .gdpar_geom_pack_fit(cs_fit, elapsed)
}

.gdpar_geom_sample_gdpar <- function(norm, n_knob, controls, seed, ...) {
  tgt <- norm$gdpar_target
  idx <- sample.int(norm$n_full, size = min(n_knob, norm$n_full))
  sub <- tgt$data[idx, , drop = FALSE]
  args <- tgt
  args$type <- NULL
  args$data <- sub
  args$iter_warmup <- controls$warmup
  args$iter_sampling <- controls$sampling
  args$chains <- controls$chains
  args$adapt_delta <- controls$adapt_delta
  args$max_treedepth <- controls$max_treedepth
  args$seed <- seed
  args$refresh <- 0L
  args$verbose <- FALSE
  args$skip_id_check <- TRUE
  t0 <- Sys.time()
  fit <- do.call(gdpar, c(args, list(...)))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .gdpar_geom_pack_fit(fit$fit, elapsed)
}

.gdpar_geom_pack_fit <- function(cs_fit, elapsed) {
  draws_arr <- cs_fit$draws()
  vars <- posterior::variables(draws_arr)
  param_vars <- setdiff(vars, "lp__")
  pm <- posterior::as_draws_matrix(
    posterior::subset_draws(draws_arr, variable = param_vars))
  pm <- as.matrix(pm)
  chain_id <- posterior::as_draws_df(draws_arr)$.chain
  list(cs_fit = cs_fit, pm = pm, chain_id = chain_id, elapsed = elapsed)
}

# Compile a Stan program once; cmdstanr caches by file hash, so repeated
# diagnostic calls and pilots at different sizes reuse the compiled model.
.gdpar_geom_compile <- function(stan_code) {
  f <- cmdstanr::write_stan_file(stan_code)
  cmdstanr::cmdstan_model(f)
}

# --- size-invariant signal extraction --------------------------------------
.gdpar_geom_extract_signals <- function(cs_fit, pm, chain_id, bounds,
                                        max_treedepth, elapsed) {
  sd_df <- cs_fit$sampler_diagnostics(format = "draws_df")
  ds <- cs_fit$diagnostic_summary(quiet = TRUE)
  total <- nrow(sd_df)
  divergent_rate <- sum(sd_df$divergent__) / total
  treedepth_sat_rate <- mean(sd_df$treedepth__ >= max_treedepth)
  ebfmi_min <- suppressWarnings(min(ds$ebfmi))
  mean_leapfrog <- mean(sd_df$n_leapfrog__)
  step_size <- mean(sd_df$stepsize__)

  cv <- stats::cov(pm)
  ev <- tryCatch(eigen(cv, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  lambda_max <- if (all(is.na(ev))) NA_real_ else max(ev)
  lambda_min <- if (all(is.na(ev))) NA_real_ else max(min(ev), 0)
  condition_number <- if (is.na(lambda_max)) NA_real_ else
    lambda_max / max(lambda_min, .Machine$double.eps)
  marg_sd <- sqrt(pmax(diag(cv), .Machine$double.eps))
  step_scale_ratio <- step_size / min(marg_sd)

  multimodality <- .gdpar_geom_multimodality(pm, chain_id)
  heavy_kurtosis <- .gdpar_geom_heavy_kurtosis(pm)
  boundary_proximity <- .gdpar_geom_boundary_prox(cs_fit, bounds)

  list(
    divergent_rate = divergent_rate, ebfmi_min = ebfmi_min,
    treedepth_sat_rate = treedepth_sat_rate,
    condition_number = condition_number, lambda_max_cov = lambda_max,
    step_scale_ratio = step_scale_ratio, mean_leapfrog = mean_leapfrog,
    multimodality = multimodality, heavy_kurtosis = heavy_kurtosis,
    boundary_proximity = boundary_proximity, elapsed = elapsed,
    n_sampling = total
  )
}

.gdpar_geom_na_signals <- function() {
  list(divergent_rate = NA_real_, ebfmi_min = NA_real_,
       treedepth_sat_rate = NA_real_, condition_number = NA_real_,
       lambda_max_cov = NA_real_, step_scale_ratio = NA_real_,
       mean_leapfrog = NA_real_, multimodality = NA_real_,
       heavy_kurtosis = NA_real_, boundary_proximity = NA_real_,
       elapsed = NA_real_, n_sampling = NA_real_)
}

# Between-chain mean separation relative to within-chain spread: chains that
# land in different modes inflate this. Robust to a short pilot.
.gdpar_geom_multimodality <- function(pm, chain_id) {
  chains <- unique(chain_id)
  if (length(chains) < 2L) return(0)
  scores <- vapply(seq_len(ncol(pm)), function(j) {
    cmeans <- vapply(chains, function(ch) mean(pm[chain_id == ch, j]),
                     numeric(1))
    csds <- vapply(chains, function(ch) stats::sd(pm[chain_id == ch, j]),
                   numeric(1))
    within <- mean(csds, na.rm = TRUE)
    if (!is.finite(within) || within <= 0) return(0)
    stats::sd(cmeans) / within
  }, numeric(1))
  max(scores, na.rm = TRUE)
}

# Maximum sample excess kurtosis across coordinates (heavy tails).
.gdpar_geom_heavy_kurtosis <- function(pm) {
  ks <- vapply(seq_len(ncol(pm)), function(j) {
    x <- pm[, j]
    m <- mean(x)
    s2 <- mean((x - m)^2)
    if (s2 <= 0) return(0)
    mean((x - m)^4) / s2^2 - 3
  }, numeric(1))
  max(ks, na.rm = TRUE)
}

# Fraction of constrained draws within boundary_eps (relative to the bound
# range) of a declared bound, maximised over bounded parameters.
.gdpar_geom_boundary_prox <- function(cs_fit, bounds, eps = 0.01) {
  if (is.null(bounds) || length(bounds) == 0L) return(0)
  draws_arr <- cs_fit$draws()
  vars <- posterior::variables(draws_arr)
  out <- vapply(names(bounds), function(nm) {
    if (!nm %in% vars) return(0)
    x <- as.numeric(posterior::as_draws_matrix(
      posterior::subset_draws(draws_arr, variable = nm)))
    lo <- bounds[[nm]][1]
    hi <- bounds[[nm]][2]
    rng <- hi - lo
    near_lo <- mean(x <= lo + eps * rng)
    near_hi <- mean(x >= hi - eps * rng)
    max(near_lo, near_hi)
  }, numeric(1))
  max(out, na.rm = TRUE)
}

# --- culprit localisation --------------------------------------------------
.gdpar_geom_empty_culprit <- function() {
  data.frame(parameter = character(0), mechanism = character(0),
             score = numeric(0), stringsAsFactors = FALSE)
}

.gdpar_geom_localize <- function(cs_fit, pm, chain_id, bounds) {
  rows <- list()
  pnames <- colnames(pm)
  # (a) flat direction: the largest-variance eigenvector loadings.
  cv <- stats::cov(pm)
  eg <- tryCatch(eigen(cv, symmetric = TRUE), error = function(e) NULL)
  if (!is.null(eg)) {
    load <- abs(eg$vectors[, 1])
    ord <- order(load, decreasing = TRUE)
    top <- ord[seq_len(min(2L, length(ord)))]
    for (k in top) {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = pnames[k], mechanism = "flat_or_anisotropic_direction",
        score = load[k] * eg$values[1], stringsAsFactors = FALSE)
    }
  }
  # (b) funnel neck: parameter whose value most separates divergent transitions.
  sd_df <- cs_fit$sampler_diagnostics(format = "draws_df")
  div <- as.logical(sd_df$divergent__)
  if (any(div) && !all(div)) {
    for (j in seq_len(ncol(pm))) {
      x <- pm[, j]
      sdx <- stats::sd(x)
      if (!is.finite(sdx) || sdx <= 0) next
      shift <- abs(mean(x[div]) - mean(x[!div])) / sdx
      if (shift > 0.5) {
        rows[[length(rows) + 1L]] <- data.frame(
          parameter = pnames[j], mechanism = "divergence_neck",
          score = shift, stringsAsFactors = FALSE)
      }
    }
  }
  # (c) boundary: bounded parameters piling at an edge.
  if (!is.null(bounds)) {
    for (nm in names(bounds)) {
      prox <- .gdpar_geom_boundary_prox(cs_fit, bounds[nm])
      if (prox > 0.05) {
        rows[[length(rows) + 1L]] <- data.frame(
          parameter = nm, mechanism = "boundary_pile",
          score = prox, stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows) == 0L) return(.gdpar_geom_empty_culprit())
  out <- do.call(rbind, rows)
  out <- out[order(out$score, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# --- difficulty-vs-n curve -------------------------------------------------
.gdpar_geom_difficulty_curve <- function(signals, thresholds) {
  ok <- !signals$failed & is.finite(signals$condition_number) &
    is.finite(signals$n)
  s <- signals[ok, , drop = FALSE]
  if (nrow(s) < 2L || length(unique(s$n)) < 2L) {
    return(list(slope = NA_real_, grows_with_n = NA,
                note = "need at least two distinct sizes"))
  }
  fit <- stats::lm(log(condition_number) ~ log(n), data = s)
  slope <- unname(stats::coef(fit)[2])
  list(slope = slope,
       grows_with_n = slope > thresholds$nslope_grows,
       note = NULL)
}

# --- rule-based classifier -------------------------------------------------
.gdpar_geom_remedy_for <- function(pathology) {
  switch(pathology,
    isotropic           = list(remedy = "euclidean_diagonal", level = 0L),
    anisotropic         = list(remedy = "euclidean_dense", level = 1L),
    funnel              = list(remedy = "riemannian", level = 3L),
    heavy_tails         = list(remedy = "finsler_relativistic", level = 4L),
    quasi_deterministic = list(remedy = "sub_riemannian", level = 5L),
    multimodal          = list(remedy = "tempering", level = 6L),
    boundary            = list(remedy = "boundary_reparam", level = 6L),
    flat_direction      = list(remedy = "reparam_eliminate", level = -1L),
    list(remedy = "unknown", level = NA_integer_))
}

.gdpar_geom_classify <- function(sig, ncurve, culprit, th) {
  trace <- character(0)
  decide <- function(p, conf) {
    rem <- .gdpar_geom_remedy_for(p)
    list(pathology = p, confidence = conf, remedy = rem$remedy,
         level = rem$level, trace = trace)
  }
  grows <- isTRUE(ncurve$grows_with_n)
  sq <- function(x, ref) {
    if (!is.finite(x) || !is.finite(ref) || ref <= 0) return(0.5)
    min(1, max(0, x / ref / 2 + 0.25))
  }
  kurt_high <- isTRUE(sig$heavy_kurtosis >= th$heavy_kurtosis_high)
  cond_high <- isTRUE(sig$condition_number >= th$condition_high)
  # 1. boundary: bounded parameters piling at an edge.
  if (isTRUE(sig$boundary_proximity >= th$boundary_prox_high)) {
    trace <- c(trace, "boundary_proximity >= threshold")
    return(decide("boundary", sq(sig$boundary_proximity, th$boundary_prox_high)))
  }
  # 2. multimodal: chains separate into distinct modes.
  if (isTRUE(sig$multimodality >= th$multimodal_high)) {
    trace <- c(trace, "multimodality >= threshold")
    return(decide("multimodal", sq(sig$multimodality, th$multimodal_high)))
  }
  # 3. flat direction: a near-zero Hessian eigenvalue, i.e. a posterior
  #    direction with prior-scale (huge) variance that does not shrink with n.
  if (isTRUE(sig$lambda_max_cov >= th$flat_var_high) && !grows) {
    trace <- c(trace, "lambda_max_cov >= threshold and not n-growing")
    return(decide("flat_direction", sq(sig$lambda_max_cov, th$flat_var_high)))
  }
  # 4. quasi-deterministic: high condition number that grows with n
  #    (smallest direction pinned ever tighter; the eBird count case).
  if (cond_high && grows) {
    trace <- c(trace, "condition_number high and grows with n")
    return(decide("quasi_deterministic",
                  sq(sig$condition_number, th$condition_high)))
  }
  # 5. anisotropic: high but n-flat condition number, energy and tails benign.
  if (cond_high && !kurt_high) {
    trace <- c(trace, "condition_number high, n-flat, benign tails")
    return(decide("anisotropic", sq(sig$condition_number, th$condition_high)))
  }
  # 6. funnel: variable curvature -- heavy marginals (a neck-driven kurtosis)
  #    together with a strongly depressed energy fraction, OR an explicit
  #    divergence neck with tiny steps.
  funnel_energy <- kurt_high && isTRUE(sig$ebfmi_min <= th$funnel_ebfmi_low)
  funnel_div <- isTRUE(sig$divergent_rate >= th$divergent_rate_high) &&
    isTRUE(sig$step_scale_ratio <= th$step_scale_ratio_low)
  if (funnel_energy || funnel_div) {
    trace <- c(trace,
               if (funnel_energy) "high kurtosis with very low E-BFMI"
               else "divergence neck with tiny step/scale ratio")
    return(decide("funnel", if (funnel_energy)
      sq(sig$heavy_kurtosis, th$heavy_kurtosis_high) else
        sq(sig$divergent_rate, th$divergent_rate_high)))
  }
  # 7. heavy tails: heavy marginals with an isotropic (low-condition) geometry.
  if (kurt_high && isTRUE(sig$condition_number <= th$heavy_cond_max)) {
    trace <- c(trace, "high kurtosis with low condition number")
    return(decide("heavy_tails", sq(sig$heavy_kurtosis, th$heavy_kurtosis_high)))
  }
  # 8. residual high condition (n-curve unknown): treat as anisotropic.
  if (cond_high) {
    trace <- c(trace, "condition_number high, n-curve unknown")
    return(decide("anisotropic", 0.4))
  }
  # 9. isotropic (benign): the default Euclidean metric suffices.
  trace <- c(trace, "no pathology threshold exceeded")
  decide("isotropic", 0.6)
}

# --- cost estimation (E-BFMI / ESS-free) -----------------------------------
.gdpar_geom_cost <- function(pilot_max, ncurve) {
  sig <- pilot_max$signals
  if (isTRUE(pilot_max$failed) || !is.finite(sig$elapsed)) {
    return(list(seconds_per_1000_draws = NA_real_, mean_leapfrog = NA_real_,
                treedepth_saturation = NA_real_, tractability = "unknown"))
  }
  per_1000 <- (sig$elapsed / max(sig$n_sampling, 1)) * 1000
  sat <- sig$treedepth_sat_rate
  grows <- isTRUE(ncurve$grows_with_n)
  tract <- if (isTRUE(sat >= 0.5) && grows) {
    "intractable (escalate geometry or certify limit)"
  } else if (isTRUE(sat >= 0.5)) {
    "expensive"
  } else {
    "tractable"
  }
  list(seconds_per_1000_draws = per_1000, mean_leapfrog = sig$mean_leapfrog,
       treedepth_saturation = sat, tractability = tract)
}

#' Print method for gdpar_geometry_diagnostic objects
#'
#' @param x An object of class \code{gdpar_geometry_diagnostic}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geometry_diagnostic <- function(x, ...) {
  cat("<gdpar_geometry_diagnostic>\n")
  cat("  pathology:   ", x$pathology,
      "  (confidence ", format(x$confidence, digits = 2), ")\n", sep = "")
  cat("  remedy:      ", x$recommended_geometry,
      "  (geometry level ", x$geometry_level, ")\n", sep = "")
  if (!is.null(x$difficulty_curve$slope)) {
    cat("  n-curve:     slope(log kappa ~ log n) = ",
        format(x$difficulty_curve$slope, digits = 3),
        " | grows with n: ", x$difficulty_curve$grows_with_n, "\n", sep = "")
  }
  cat("  cost:        ", format(x$cost$seconds_per_1000_draws, digits = 3),
      " s / 1000 draws | tractability: ", x$cost$tractability, "\n", sep = "")
  if (nrow(x$culprit) > 0L) {
    cat("  culprit(s):  ",
        paste(utils::head(x$culprit$parameter, 3), collapse = ", "), "\n",
        sep = "")
  }
  if (!is.null(x$correct)) {
    cat("  ground truth:", x$ground_truth$pathology,
        " | classified correctly:", x$correct, "\n")
  }
  invisible(x)
}

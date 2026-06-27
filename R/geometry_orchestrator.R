# Geometry-adaptive sampling orchestrator (Block RG, RG.5).
#
# This file closes the loop of the Block RG charter (section 2.4): an opt-in
# controller that DIAGNOSES the geometry of a posterior (reusing
# gdpar_geometry_diagnostic), CLASSIFIES the pathology, SELECTS a level of the
# sampler hierarchy {euclidean diagonal -> dense -> Riemannian -> relativistic
# -> sub-Riemannian}, SAMPLES with that level (reusing gdpar_geom_hmc and the
# RG.2/RG.3/RG.4 metrics), RE-DIAGNOSES, and either ESCALATES the level or emits
# a CERTIFIED LIMIT -- a first-class, reproducible, falsifiable output recording
# why every attempted level failed within the budget, plus a conjectured
# prescription of the smallest change predicted to break the limit.
#
# Design decisions (session B9.29, canonised D95). Four hinges were consulted
# upfront and decided by the user under the cornerstone rule (maximum
# multidimensional robustness, cost irrelevant):
#
#   (a) Level selection -- the COMBINED map. The diagnostic's transparent,
#       calibrated rule-based classifier is primary; a continuous proximity
#       score over the size-invariant signals is computed as a SECOND, independent
#       estimator. When the two agree and confidence is high, the discrete level
#       is used directly (full transparency). When they disagree, or confidence
#       is low, the controller starts conservatively at the LOWER of the two
#       candidate ladder levels and lets the escalation climb -- robust precisely
#       at the funnel/heavy-tail border that RG.1.c found mutually confusable
#       (recall 0.63/0.71). A user override (level_map / entry_level) sits on top.
#
#   (b) Escalation -- the closed-loop, RE-DIAGNOSIS-guided controller, ARMOURED.
#       On failure the controller re-diagnoses with a fresh pilot and lets the
#       NEW signature pick the next level (a re-diagnosis-guided upward jump that
#       may skip ladder rungs), under a
#       monotone ratchet, a visited-state memo (provably acyclic walk), per-level
#       caps, a global round cap, a no-progress (stall) detector, hierarchical
#       budget accounting with cost-aware admission and cheap-probe-then-full
#       tiering (successive-halving spirit), a per-fit wall-time watchdog,
#       graceful degradation (always return the best result so far), deterministic
#       seeding (the whole adaptive trajectory is bit-reproducible despite being
#       data-adaptive, so a re-run retraces it identically -- that determinism is
#       the resumability guarantee, with an optional atomic ledger checkpoint for
#       durability), a multi-signal success gate with hysteresis, a frozen-level
#       final sampling phase (preserving exactness/ergodicity of the returned
#       draws), and the capability-subsumption guarantee (level L+1 strictly
#       generalises level L, so the monotone ladder fallback is a proven
#       monotone-improvement backstop when an adaptive jump misfires). Every level
#       is Metropolis-exact, so a wrong jump never corrupts correctness, only
#       wastes a bounded amount of budget.
#
#   (c) Certified limit -- the EVIDENCE LEDGER plus a CONJECTURED PRESCRIPTION.
#       The certificate separates, in three rigour layers (algebraic / statistical
#       / numerical, never conflated), the DEMONSTRATED record (what was measured)
#       from the conjectured prescription (the smallest change the evidence
#       predicts would succeed), the latter tagged as a conjecture in the sense of
#       ORPHEUS-PIMC section 16.3 and made falsifiable (a re-run test). It mirrors
#       the diagnostic's reproducibility block.
#
#   (d) Scope -- STANDALONE with a rich, structured return. The controller never
#       touches gdpar()'s fit path (the default branch is bit-identical, goldens
#       intact); the return carries everything RG.6 will need to wire it into a
#       real fit (the winning metric, the ledger, the reproducibility block)
#       WITHOUT freezing a contract that RG.6 might have to redesign. The
#       statistical-layer validation exercises a real compiled cmdstan model end
#       to end (the B9.26/B9.27 bridge).
#
# Anchors (gdpar integrates established patterns, it does not claim to invent
# them): Betancourt (E-BFMI, sampling pathologies); Girolami & Calderhead 2011;
# Roberts & Rosenthal 2007 (diminishing adaptation / containment, the stability
# rationale of the level sequence); Li et al. 2018 (Hyperband / successive
# halving, the cheap-probe-then-full budget tiering); Xu et al. 2008 (SATzilla,
# runtime-predicted algorithm selection with a presolver fallback -- here the
# monotone ladder).

# ---------------------------------------------------------------------------
# Tunable data (exposed, not hard-coded, so they can be re-calibrated).
# ---------------------------------------------------------------------------

#' Success criteria for the geometry-adaptive orchestrator
#'
#' The multi-signal success gate of \code{\link{gdpar_geom_orchestrate}}: a level
#' is accepted only when a \emph{conjunction} of size-invariant sampler signals
#' holds, never a single number. Using a conjunction (and never R-hat or the
#' effective sample size on short runs) is the generalisation of the
#' \code{rhat = Inf} false-positive lesson of sessions B9.20/B9.21: a sampler
#' that accepts every proposal but never moves must not be mistaken for success.
#'
#' @details
#' The gate requires the acceptance rate inside a healthy band, the divergence
#' rate below a ceiling, and the energy Bayesian fraction of missing information
#' (\dQuote{E-BFMI}) above a floor. The decisive (full-budget) gate is applied
#' with a hysteresis margin (the thresholds tightened by \code{1 + margin}) so a
#' level must \emph{clearly} pass to be declared a success, preventing chattering
#' around the boundary.
#'
#' @return A named list of numeric criteria: \code{accept_low}, \code{accept_high}
#'   (the healthy acceptance band), \code{divergent_rate_high} (the divergence
#'   ceiling) and \code{ebfmi_low} (the E-BFMI floor).
#' @seealso \code{\link{gdpar_geom_orchestrate}},
#'   \code{\link{gdpar_geom_orchestrate_budget}}.
#' @examples
#' str(gdpar_geom_orchestrate_criteria())
#' @export
gdpar_geom_orchestrate_criteria <- function() {
  list(
    accept_low          = 0.5,
    accept_high         = 0.999,
    divergent_rate_high = 0.02,
    ebfmi_low           = 0.3
  )
}

#' Budget and stopping rule for the geometry-adaptive orchestrator
#'
#' The budget that bounds \code{\link{gdpar_geom_orchestrate}}'s closed loop and
#' the knobs of its armour. The loop spends fits and wall-clock time against
#' these caps, checks the remaining budget \emph{before} starting each fit (no
#' half-fits), and emits a certified limit rather than overrunning.
#'
#' @details
#' The fields are: \code{max_rounds} (hard ceiling on diagnose--sample cycles),
#' \code{max_levels} (cap on distinct ladder levels tried), \code{probe_warmup} /
#' \code{probe_iter} (the cheap probe budget every level gets first) and
#' \code{full_warmup} / \code{full_iter} (the full budget a level earns only by
#' passing its probe -- the successive-halving spirit), \code{epsilon} / \code{L}
#' (the base leapfrog step and the trajectory length), \code{tune_epsilon}
#' (whether to run a short coarse step-size search per level before the probe,
#' targeting a healthy acceptance band, so a level is not failed merely because
#' the base step was ill-matched to its geometry; the sub-Riemannian level
#' searches around its own \code{suggested_epsilon}) and \code{tune_iter}
#' (iterations per tuning probe), \code{max_seconds}
#' (global wall-time cap), \code{max_seconds_per_fit} (the per-fit watchdog),
#' \code{max_fits} (global fit cap), \code{n_rediagnose} (pilots per
#' re-diagnosis; \eqn{\ge 2} enables the stability check that falls back to the
#' safe ladder step when the re-diagnosis itself is unstable), \code{stall_limit}
#' (consecutive no-progress rounds before stopping) and \code{hysteresis} (the
#' margin tightening the decisive success gate).
#'
#' @return A named list of budget and stopping-rule settings.
#' @seealso \code{\link{gdpar_geom_orchestrate}},
#'   \code{\link{gdpar_geom_orchestrate_criteria}}.
#' @examples
#' str(gdpar_geom_orchestrate_budget())
#' @export
gdpar_geom_orchestrate_budget <- function() {
  list(
    max_rounds          = 8L,
    max_levels          = 5L,
    probe_warmup        = 150L,
    probe_iter          = 150L,
    full_warmup         = 500L,
    full_iter           = 500L,
    epsilon             = 0.25,
    L                   = 25L,
    tune_epsilon        = TRUE,
    tune_iter           = 60L,
    max_seconds         = Inf,
    max_seconds_per_fit = Inf,
    max_fits            = 40L,
    n_rediagnose        = 1L,
    stall_limit         = 2L,
    hysteresis          = 0.1
  )
}

# ---------------------------------------------------------------------------
# The sampler hierarchy (the ladder) and the pathology -> level map.
# ---------------------------------------------------------------------------

# Ordered ladder of the sampler levels the orchestrator can build, by hierarchy
# index. Levels 6 (tempering / boundary reparametrisation) and -1
# (reparametrise / eliminate, the Option A case) are deliberately OUT OF the
# orchestrator's scope: they are not sampler levels of the geometry hierarchy
# but separate remedies, so a diagnosis pointing at them short-circuits to a
# certificate whose prescription names the proper remedy (no overreach).
.gdpar_orch_ladder <- function() {
  list(
    euclidean_diagonal = list(key = "euclidean_diagonal", index = 0L,
                              needs_fisher = FALSE),
    euclidean_dense    = list(key = "euclidean_dense", index = 1L,
                              needs_fisher = FALSE),
    riemannian         = list(key = "riemannian", index = 3L,
                              needs_fisher = FALSE),
    relativistic       = list(key = "relativistic", index = 4L,
                              needs_fisher = FALSE),
    sub_riemannian     = list(key = "sub_riemannian", index = 5L,
                              needs_fisher = TRUE)
  )
}

# The in-ladder ladder keys in ascending hierarchy order.
.gdpar_orch_ladder_keys <- function() {
  c("euclidean_diagonal", "euclidean_dense", "riemannian", "relativistic",
    "sub_riemannian")
}

# Map a classified pathology to a ladder key, or NA for an out-of-scope remedy.
.gdpar_orch_pathology_to_key <- function(pathology) {
  switch(pathology,
    isotropic           = "euclidean_diagonal",
    anisotropic         = "euclidean_dense",
    funnel              = "riemannian",
    heavy_tails         = "relativistic",
    quasi_deterministic = "sub_riemannian",
    NA_character_)   # multimodal / boundary / flat_direction / unknown.
}

# The first ladder key strictly above a hierarchy index, or NA if at the top.
.gdpar_orch_next_above <- function(index) {
  ladder <- .gdpar_orch_ladder()
  keys <- .gdpar_orch_ladder_keys()
  idx <- vapply(ladder[keys], function(l) l$index, integer(1))
  above <- keys[idx > index]
  if (length(above) == 0L) NA_character_ else above[which.min(idx[idx > index])]
}

# ---------------------------------------------------------------------------
# (a) Combined level selection: discrete classifier + proximity cross-check.
# ---------------------------------------------------------------------------

# A transparent, threshold-relative proximity score per in-ladder pathology --
# the SECOND, independent estimator of the combined map. Each score is built
# from the same size-invariant signals the rule-based classifier uses, scaled by
# the classifier's own thresholds (so "proximity in standardised signal space"
# means distance relative to the calibrated cut points). The argmax is the
# proximity pathology, cross-checked against the discrete classifier's pick.
.gdpar_orch_pathology_scores <- function(sig, grows, th) {
  safe <- function(x, default = 0) if (is.finite(x)) x else default
  cond  <- safe(sig$condition_number, 1)
  kurt  <- safe(sig$heavy_kurtosis, 0)
  ebfmi <- safe(sig$ebfmi_min, 1)
  r_cond <- cond / th$condition_high
  r_kurt <- kurt / th$heavy_kurtosis_high
  ebfmi_pos <- max(ebfmi, 1e-3)
  # Funnel and heavy tails share a heavy-kurtosis signature; the calibrated
  # classifier splits them by the energy fraction (a depressed E-BFMI is the
  # funnel signature) and the condition number (a low condition with healthy
  # energy is heavy tails). The proximity scores mirror that split: the heavy-
  # tail score is damped when the E-BFMI is depressed, so a deep funnel is not
  # mistaken for heavy tails.
  c(
    isotropic           = -max(r_cond, r_kurt),
    anisotropic         = r_cond * (if (isTRUE(grows)) 0 else 1) *
                          (if (r_kurt < 1) 1 else 0.25),
    funnel              = r_kurt * (th$funnel_ebfmi_low / ebfmi_pos),
    heavy_tails         = r_kurt * (th$heavy_cond_max / max(cond, 1)) *
                          min(1, ebfmi_pos / th$funnel_ebfmi_low),
    quasi_deterministic = r_cond * (if (isTRUE(grows)) 2 else 0.5)
  )
}

# Resolve the entry level from a diagnostic (the COMBINED map). Returns the entry
# ladder key (or NA if out of scope) plus the audit fields (the discrete
# pathology, the proximity pathology, the conflict flag) carried into the ledger
# and certificate. `diag` is a gdpar_geometry_diagnostic (or a compatible list).
.gdpar_orch_select_entry <- function(diag, level_map, entry_level, th) {
  ladder <- .gdpar_orch_ladder()
  if (!is.null(entry_level)) {
    if (!entry_level %in% names(ladder)) {
      gdpar_abort(
        sprintf("'entry_level' must be one of: %s.",
                paste(names(ladder), collapse = ", ")),
        class = "gdpar_input_error")
    }
    return(list(key = entry_level, discrete = diag$pathology,
                proximity = NA_character_, conflict = FALSE,
                source = "user_entry_level"))
  }
  discrete <- diag$pathology
  if (!is.null(level_map) && !is.null(level_map[[discrete]])) {
    return(list(key = level_map[[discrete]], discrete = discrete,
                proximity = NA_character_, conflict = FALSE,
                source = "user_level_map"))
  }
  discrete_key <- .gdpar_orch_pathology_to_key(discrete)
  # The proximity second estimator over the size-invariant signals.
  sig <- diag$signals[nrow(diag$signals), , drop = FALSE]
  grows <- isTRUE(diag$difficulty_curve$grows_with_n)
  scores <- .gdpar_orch_pathology_scores(sig, grows, th)
  proximity <- names(scores)[which.max(scores)]
  proximity_key <- .gdpar_orch_pathology_to_key(proximity)
  conflict <- !identical(proximity, discrete)
  confident <- isTRUE(diag$confidence >= 0.6)
  # If discrete points out of scope, defer to that (handled by the caller as a
  # short-circuit certificate). Otherwise, on a conflict or low confidence, start
  # conservatively at the lower of the two candidate ladder levels and let the
  # escalation climb; else use the discrete level directly.
  if (is.na(discrete_key)) {
    key <- NA_character_
  } else if ((conflict || !confident) && !is.na(proximity_key)) {
    idx_d <- ladder[[discrete_key]]$index
    idx_p <- ladder[[proximity_key]]$index
    key <- if (idx_p < idx_d) proximity_key else discrete_key
  } else {
    key <- discrete_key
  }
  list(key = key, discrete = discrete, proximity = proximity,
       conflict = conflict, confident = confident, source = "combined_map")
}

# ---------------------------------------------------------------------------
# Metric construction per level (with per-level guards).
# ---------------------------------------------------------------------------

# A constant dense Euclidean preconditioner: the expected Fisher at the reference
# when supplied, else the precision from warmup draws, else the SoftAbs curvature
# at the reference (always available). Floored to be strictly positive-definite.
.gdpar_orch_dense_mass <- function(geom_target, fisher, reference, warmup_draws) {
  d <- geom_target$dim
  M <- NULL
  if (is.function(fisher)) {
    M <- tryCatch(as.matrix(fisher(reference)), error = function(e) NULL)
  }
  if (is.null(M) && is.matrix(warmup_draws) && nrow(warmup_draws) > d) {
    cv <- tryCatch(stats::cov(warmup_draws), error = function(e) NULL)
    if (!is.null(cv)) {
      M <- tryCatch(chol2inv(.gdpar_geom_chol_spd(cv)), error = function(e) NULL)
    }
  }
  if (is.null(M)) {
    soft <- gdpar_geom_metric_riemannian(geom_target, curvature = "softabs")
    M <- soft$mass(reference)
  }
  .gdpar_geom_floor_spd(as.matrix(M), 1e-8)
}

# Build the metric for a ladder key; a malformed metric is reported (not a crash)
# so the controller can record it and escalate.
.gdpar_orch_build_metric <- function(key, geom_target, fisher, reference,
                                     speed, rest_mass, warmup_draws) {
  d <- geom_target$dim
  out <- tryCatch({
    metric <- switch(key,
      euclidean_diagonal = gdpar_geom_metric_euclidean(dim = d),
      euclidean_dense = gdpar_geom_metric_euclidean(
        M = .gdpar_orch_dense_mass(geom_target, fisher, reference, warmup_draws)),
      riemannian = gdpar_geom_metric_riemannian(
        geom_target,
        curvature = if (is.function(fisher)) "fisher" else "softabs",
        fisher = fisher),
      relativistic = gdpar_geom_metric_relativistic(
        geom_target,
        curvature = if (is.function(fisher)) "fisher" else "softabs",
        fisher = fisher, speed = speed, rest_mass = rest_mass),
      sub_riemannian = {
        if (!is.function(fisher)) {
          stop("the sub-Riemannian level requires an expected-Fisher function")
        }
        gdpar_geom_metric_subriemannian(geom_target, fisher = fisher,
                                        reference = reference)
      },
      stop(sprintf("unknown ladder key '%s'", key)))
    list(ok = TRUE, metric = metric, reason = NA_character_)
  }, error = function(e) list(ok = FALSE, metric = NULL,
                              reason = conditionMessage(e)))
  out
}

# ---------------------------------------------------------------------------
# The multi-signal success gate (with hysteresis).
# ---------------------------------------------------------------------------

# Evaluate the gate on one fit. `total` is the number of proposals (warmup +
# sampling) the divergence rate is taken against. `margin` tightens the
# thresholds (hysteresis) for the decisive full gate; pass 0 for the probe gate.
.gdpar_orch_success_gate <- function(fit, total, criteria, margin = 0) {
  if (is.null(fit) || is.null(fit$accept_rate)) {
    return(list(pass = FALSE, accept = NA_real_, divergent_rate = NA_real_,
                ebfmi = NA_real_, reasons = "no fit"))
  }
  m <- 1 + margin
  accept <- fit$accept_rate
  div_rate <- fit$n_divergent / max(total, 1)
  ebfmi <- fit$ebfmi
  pass_accept <- isTRUE(accept >= criteria$accept_low * m) &&
    isTRUE(accept <= criteria$accept_high)
  pass_div <- isTRUE(div_rate <= criteria$divergent_rate_high / m)
  pass_ebfmi <- is.finite(ebfmi) && isTRUE(ebfmi >= criteria$ebfmi_low * m)
  reasons <- c(
    if (!pass_accept) "acceptance outside healthy band",
    if (!pass_div) "divergence rate too high",
    if (!pass_ebfmi) "E-BFMI below floor")
  list(pass = pass_accept && pass_div && pass_ebfmi,
       accept = accept, divergent_rate = div_rate, ebfmi = ebfmi,
       reasons = if (length(reasons)) reasons else "all criteria met")
}

# A scalar "health" score for best-so-far tracking and the no-progress detector;
# higher is better. Rewards E-BFMI and acceptance, penalises divergences.
.gdpar_orch_score <- function(fit, total) {
  if (is.null(fit) || is.null(fit$accept_rate)) return(-Inf)
  ebfmi <- if (is.finite(fit$ebfmi)) fit$ebfmi else 0
  div_rate <- fit$n_divergent / max(total, 1)
  ebfmi + 0.2 * fit$accept_rate - 5 * div_rate
}

# A short coarse step-size search: the orchestrator selects the GEOMETRY, but a
# level must not be failed merely because the base step was ill-matched to it.
# A handful of cheap probes over a geometric epsilon grid (around the level's
# own suggested step where it has one) pick the step whose acceptance is closest
# to a healthy target, penalising divergences -- a deterministic, bounded
# analogue of step-size adaptation. Returns the chosen step and the fits spent.
.gdpar_orch_tune_epsilon <- function(geom_target, metric, budget, base_seed,
                                     round) {
  base <- if (!is.null(metric$suggested_epsilon)) metric$suggested_epsilon
          else budget$epsilon
  grid <- sort(unique(pmin(base * c(2, 1, 0.5, 0.25, 0.1), 5)))
  target_acc <- 0.8
  it <- budget$tune_iter
  best_eps <- base; best_gap <- Inf; nfit <- 0L; secs <- 0
  for (k in seq_along(grid)) {
    eps <- grid[k]
    t0 <- Sys.time()
    fit <- tryCatch(
      gdpar_geom_hmc(geom_target, metric = metric, epsilon = eps, L = budget$L,
                     n_iter = it, n_warmup = it,
                     seed = .gdpar_orch_seed(base_seed, round, 100L + k)),
      error = function(e) NULL)
    secs <- secs + as.numeric(difftime(Sys.time(), t0, units = "secs"))
    nfit <- nfit + 1L
    if (is.null(fit)) next
    div_rate <- fit$n_divergent / (2 * it)
    gap <- abs(fit$accept_rate - target_acc) + 5 * div_rate
    if (is.finite(gap) && gap < best_gap) { best_gap <- gap; best_eps <- eps }
  }
  list(epsilon = best_eps, fits = nfit, seconds = secs)
}

# ---------------------------------------------------------------------------
# Deterministic seeding: the whole adaptive trajectory is a function of the base
# seed, the round and the slot, so a re-run retraces it identically.
# ---------------------------------------------------------------------------
.gdpar_orch_seed <- function(base, round, slot) {
  as.integer((base + round * 7919L + slot * 104729L) %% 2147483646L) + 1L
}

# ---------------------------------------------------------------------------
# (c) The certified limit: evidence ledger (three layers) + conjectured
# prescription.
# ---------------------------------------------------------------------------

# Derive the conjectured, falsifiable prescription from the accumulated evidence:
# the smallest change predicted to break the limit. Each item is tagged as a
# conjecture (ORPHEUS section 16.3) with a falsifiable re-run test.
.gdpar_orch_prescription <- function(final_diag, ledger, fisher_available,
                                     budget, best_score, tried_keys) {
  items <- list()
  add <- function(action, rationale, test) {
    items[[length(items) + 1L]] <<- list(
      action = action, rationale = rationale, falsifiable_test = test,
      status = "conjectured")
  }
  patho <- final_diag$pathology
  grows <- isTRUE(final_diag$difficulty_curve$grows_with_n)
  slope <- final_diag$difficulty_curve$slope
  # Out-of-scope remedies: the geometry sampler ladder does not address these, so
  # name ONLY the proper remedy (no spurious geometry advice).
  if (identical(patho, "multimodal")) {
    add("apply tempering / a multimodal sampler (out of the geometry ladder)",
        "chains separate into distinct modes; no position metric bridges them",
        "re-run with a tempered or mode-jumping sampler")
    return(items)
  }
  if (identical(patho, "boundary")) {
    add("reparametrise the boundary-pinned parameter (out of the geometry ladder)",
        "mass pinned against a bound gives singular curvature at the edge",
        "re-run on a reparametrised, unbounded scale")
    return(items)
  }
  if (identical(patho, "flat_direction")) {
    add("reparametrise or eliminate the flat direction (Option A)",
        "a near-zero Hessian eigenvalue = non-identification, not curvature",
        "re-run after removing the unidentified scale / contrast")
    return(items)
  }
  # Was the sub-Riemannian level actually sampled (built successfully), or only
  # attempted and unbuildable (no expected Fisher)?
  sub_sampled <- any(vapply(ledger, function(r)
    identical(r$level, "sub_riemannian") && isTRUE(r$metric_ok), logical(1)))
  # Quasi-determinism reachable only with the expected Fisher: prescribe it when
  # the geometry contracts (grows with n, or the class is quasi-deterministic)
  # and the sub-Riemannian level could not be sampled for want of a Fisher.
  if (!fisher_available && !sub_sampled &&
      (grows || identical(patho, "quasi_deterministic"))) {
    add("supply the expected Fisher and enable the sub-Riemannian level",
        sprintf("the geometry contracts (difficulty-vs-n slope %.2f): the sub-Riemannian remedy needs the expected Fisher, which was not supplied",
                if (is.finite(slope)) slope else NA_real_),
        "re-run with fisher = the expected Fisher (or gdpar_geom_fisher_simulator)")
  }
  # Sub-Riemannian sampled from the origin but still failing: relocate the
  # reference to the explored mode.
  if (sub_sampled) {
    add("set the sub-Riemannian reference to the warmup mode, not the origin",
        "the wall curvature is read at the reference; an off-mode reference mis-sets the floor/wall split",
        "re-run with reference = the best position found")
  }
  # Budget: extrapolate from the cost verdict.
  tract <- final_diag$cost$tractability
  if (is.finite(best_score) && best_score > -Inf &&
      identical(tract, "expensive")) {
    add("increase the per-fit budget (more iterations) and tune epsilon / L",
        "sampling was improving but expensive (tree-depth saturation)",
        "re-run with full_iter and full_warmup scaled up by a factor")
  }
  if (length(items) == 0L) {
    add("increase the budget (rounds / iterations) and revisit epsilon / L",
        "no level met the success gate within the current budget",
        "re-run with a larger gdpar_geom_orchestrate_budget()")
  }
  items
}

# Assemble the certificate object (three rigour layers + reproducibility +
# prescription).
.gdpar_orch_certificate <- function(verdict, final_diag, ledger, culprit,
                                    prescription, reproducibility, budget_spent) {
  obj <- list(
    verdict = verdict,
    # Algebraic layer: what the geometry IS.
    algebraic = list(
      pathology = final_diag$pathology,
      condition_number =
        final_diag$signals$condition_number[nrow(final_diag$signals)],
      difficulty_slope = final_diag$difficulty_curve$slope,
      grows_with_n = final_diag$difficulty_curve$grows_with_n,
      culprit = culprit),
    # Statistical layer: how sampling behaved at each level.
    statistical = ledger,
    # Numerical layer: the computational facts (budget, fits, time).
    numerical = budget_spent,
    prescription = prescription,
    reproducibility = reproducibility)
  class(obj) <- c("gdpar_geom_certificate", "list")
  obj
}

# ---------------------------------------------------------------------------
# (RG.7 step 4, D99) The opt-in Laplace fallback.
# ---------------------------------------------------------------------------

# When a run ends in a certified limit (the genuinely non-Gaussian canyon the
# geometry ladder cannot sample) and the caller opted in, attach the Laplace
# approximation and relabel the status `certified_limit_laplace` -- the sampling
# limit STILL stands, and a labelled, fidelity-diagnosed Laplace posterior is
# provided (the mgcv/REML and INLA/Laplace regime). On the out-of-scope path the
# Gaussian premise is usually violated (multimodality, a flat direction or a
# boundary all give a non-positive-definite curvature), so the fallback is
# attached there ONLY when the curvature at the mode is genuinely PD
# (status `out_of_scope_laplace`); otherwise the certificate is left untouched
# (no overreach). The Laplace is never advertised as exact; its scalar
# fit-quality label carries the honesty. Default off => bit-identical output.
.gdpar_orch_attach_laplace <- function(obj, geom_target, reference, draws, seed,
                                       inform) {
  lap <- tryCatch(
    gdpar_geom_laplace(geom_target, reference = reference, draws = draws,
                       seed = seed),
    error = function(e) {
      gdpar_warn(sprintf(paste(
        "laplace_fallback: could not build the Laplace approximation (%s);",
        "leaving the certificate unchanged."), conditionMessage(e)),
        class = "gdpar_geometry_warning")
      NULL
    })
  if (is.null(lap)) return(obj)
  base <- obj$status
  if (identical(base, "out_of_scope") && !isTRUE(lap$all_pos)) {
    inform(paste("laplace_fallback: the out-of-scope curvature is not",
                 "positive-definite; the Gaussian premise does not hold, so no",
                 "Laplace is attached."))
    return(obj)
  }
  obj$laplace <- lap
  obj$status <- paste0(base, "_laplace")
  inform(sprintf(paste("laplace_fallback: attached a '%s' Laplace approximation;",
                       "status '%s'."), lap$fit_quality_label, obj$status))
  obj
}

# ---------------------------------------------------------------------------
# The orchestrator.
# ---------------------------------------------------------------------------

#' Geometry-adaptive sampling orchestrator (opt-in)
#'
#' Diagnose the geometry of a posterior, select a level of the Block RG sampler
#' hierarchy that remedies it, sample, re-diagnose, and either escalate the level
#' or emit a \strong{certified limit}. This is the closed loop of the Block RG
#' charter: a single opt-in entry point over the diagnostic
#' (\code{\link{gdpar_geometry_diagnostic}}) and the five geometry levels
#' (\code{\link{gdpar_geom_metric_euclidean}},
#' \code{\link{gdpar_geom_metric_riemannian}},
#' \code{\link{gdpar_geom_metric_relativistic}},
#' \code{\link{gdpar_geom_metric_subriemannian}}), all sampled by
#' \code{\link{gdpar_geom_hmc}}. It is standalone and does not touch the package's
#' fit path; the default branch is bit-identical.
#'
#' @details
#' \strong{Level selection (the combined map).} The diagnostic's transparent,
#' calibrated rule-based classifier is the primary selector; a continuous
#' proximity score over the size-invariant signals is a second, independent
#' estimator. When they agree and confidence is high, the discrete level is used
#' directly; when they conflict or confidence is low, the controller starts at
#' the lower of the two candidate levels and lets the escalation climb -- robust
#' at the funnel/heavy-tail border the calibration found mutually confusable. A
#' user \code{level_map} or \code{entry_level} overrides this.
#'
#' \strong{Escalation (closed-loop, armoured).} On a failed level the controller
#' re-diagnoses with a fresh pilot and lets the new signature pick the next
#' level, under a monotone ratchet (never below a level already tried), a
#' visited-state memo (an acyclic walk), a global round cap, per-level caps, a
#' no-progress detector, hierarchical budget accounting with cost-aware admission
#' and cheap-probe-then-full tiering, a per-fit wall-time watchdog, graceful
#' degradation (the best result so far is always returned), deterministic seeding
#' (the adaptive trajectory is bit-reproducible, so a re-run retraces it -- the
#' resumability guarantee), a multi-signal success gate with hysteresis, and a
#' frozen-level final sampling phase. Each level is Metropolis-exact, so a wrong
#' jump wastes only a bounded amount of budget, never correctness; the monotone
#' ladder is the proven-generalisation backstop (level L+1 generalises level L).
#'
#' \strong{The certified limit.} When the budget is exhausted without success the
#' controller returns a \code{\link{gdpar_geom_certificate}}: the demonstrated
#' evidence in three rigour layers (algebraic = the geometry; statistical = the
#' per-level sampler diagnostics; numerical = the budget and fits) plus a
#' \emph{conjectured} prescription -- the smallest change predicted to break the
#' limit, tagged as a conjecture and made falsifiable by a re-run test. A
#' diagnosis pointing outside the geometry ladder (multimodality, a boundary, a
#' flat direction) short-circuits to such a certificate naming the proper remedy
#' (tempering, reparametrisation, Option A), without overreaching by pretending a
#' geometry metric fixes it.
#'
#' @param target The posterior the diagnostic probes: a
#'   \code{gdpar_geometry_target} from \code{\link{gdpar_geometry_suite}}, a
#'   \code{list(model = , data = , dim = )} wrapping a compiled \pkg{cmdstanr}
#'   model, a \code{list(stan_code = , stan_data = , dim = )}, or a
#'   \code{list(type = "gdpar", formula = , data = , ...)} (see
#'   \code{\link{gdpar_geometry_diagnostic}}).
#' @param geom_target The \code{\link{gdpar_geom_target}} sampled by the engine.
#'   When \code{NULL} it is derived from \code{target} for a suite target (via its
#'   \code{make()}) or a cmdstan model (via the compiled methods); supply it
#'   explicitly for other forms.
#' @param fisher Optional \code{fisher(theta)} returning the expected Fisher
#'   information at \code{theta}. The sub-Riemannian level \emph{requires} it (the
#'   level is otherwise skipped and named in the prescription); the Riemannian and
#'   relativistic levels use it when given and fall back to the SoftAbs curvature
#'   (always available) otherwise. Pass
#'   \code{\link{gdpar_geom_fisher_simulator}} for models with no closed form.
#' @param reference Optional reference position for the position-dependent levels
#'   (dense, sub-Riemannian); defaults to zeros and is warm-started to the best
#'   position found as the loop proceeds.
#' @param level_map Optional named list overriding the pathology-to-level map
#'   (names are pathology classes, values are ladder keys).
#' @param entry_level Optional ladder key forcing the entry level (one of
#'   \code{"euclidean_diagonal"}, \code{"euclidean_dense"}, \code{"riemannian"},
#'   \code{"relativistic"}, \code{"sub_riemannian"}).
#' @param budget A budget / stopping-rule list; see
#'   \code{\link{gdpar_geom_orchestrate_budget}}.
#' @param criteria A success-gate list; see
#'   \code{\link{gdpar_geom_orchestrate_criteria}}.
#' @param speed,rest_mass The relativistic level's speed and rest mass.
#' @param n_grid Optional size grid for the diagnostic pilots (forwarded).
#' @param checkpoint_dir Optional directory; when given, the running ledger is
#'   written atomically after each round for durability.
#' @param laplace_fallback Logical; when \code{TRUE} and the run ends in a
#'   certified limit (the genuinely non-Gaussian canyon the ladder cannot
#'   sample), attach a \code{\link{gdpar_geom_laplace}} approximation to the
#'   result and relabel the status \code{"certified_limit_laplace"} -- the
#'   mgcv/REML and INLA/Laplace competitor-parity regime (Block RG, RG.7). On the
#'   out-of-scope path it is attached only when the curvature at the mode is
#'   positive-definite (status \code{"out_of_scope_laplace"}). Defaults to
#'   \code{FALSE}, which leaves the output \strong{bit-identical} (no
#'   \code{$laplace}, the original status set).
#' @param laplace_draws Number of iid Laplace draws to carry on the attached
#'   approximation when \code{laplace_fallback = TRUE} (default \code{0L}: the
#'   mode plus precision Gaussian is the approximation).
#' @param seed Integer base seed; the whole adaptive trajectory is a
#'   deterministic function of it.
#' @param verbose Logical; print an opt-in progress trace. Defaults to TRUE.
#' @param ... Forwarded to \code{\link{gdpar_geometry_diagnostic}}.
#'
#' @return A list of class \code{gdpar_geom_orchestration} with \code{status}
#'   (\code{"resolved"}, \code{"certified_limit"} or \code{"out_of_scope"}; or,
#'   under \code{laplace_fallback = TRUE}, \code{"certified_limit_laplace"} /
#'   \code{"out_of_scope_laplace"} with an attached \code{$laplace}). A
#'   resolved run carries \code{draws}, the winning \code{level} and
#'   \code{metric}, the final \code{fit} and its \code{gate}; a limit or
#'   out-of-scope run carries a \code{\link{gdpar_geom_certificate}}. Every run
#'   carries the \code{ledger} (the per-round decision and sampling record), the
#'   initial \code{diagnosis}, the \code{best} result so far, \code{budget_spent}
#'   and \code{reproducibility}.
#'
#' @seealso \code{\link{gdpar_geometry_diagnostic}}, \code{\link{gdpar_geom_hmc}},
#'   \code{\link{gdpar_geom_orchestrate_budget}},
#'   \code{\link{gdpar_geom_orchestrate_criteria}}.
#'
#' @examples
#' # The budget and criteria are plain data and can be inspected / re-tuned.
#' b <- gdpar_geom_orchestrate_budget()
#' b$max_rounds
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   suite <- gdpar_geometry_suite()
#'   # A minimal run on the isotropic control: it resolves at the cheapest
#'   # (Euclidean) level. Heavier targets and the certified-limit path are shown
#'   # in the gated tests and the Block RG vignette.
#'   b$tune_epsilon <- FALSE
#'   b$probe_iter <- 60L; b$full_iter <- 80L; b$full_warmup <- 80L
#'   res <- gdpar_geom_orchestrate(suite$G0_isotropic, n_grid = 1, budget = b,
#'                                 pilot_warmup = 80L, pilot_sampling = 80L,
#'                                 verbose = FALSE)
#'   res$status
#' }
#' }
#'
#' @export
gdpar_geom_orchestrate <- function(target, geom_target = NULL, fisher = NULL,
                                   reference = NULL, level_map = NULL,
                                   entry_level = NULL, budget = NULL,
                                   criteria = NULL, speed = 10, rest_mass = 1,
                                   n_grid = NULL, checkpoint_dir = NULL,
                                   laplace_fallback = FALSE, laplace_draws = 0L,
                                   seed = 20260603L, verbose = TRUE, ...) {
  if (is.null(budget)) budget <- gdpar_geom_orchestrate_budget()
  if (is.null(criteria)) criteria <- gdpar_geom_orchestrate_criteria()
  assert_count(seed, "seed")
  assert_numeric_scalar(speed, "speed", lower = 0)
  assert_numeric_scalar(rest_mass, "rest_mass", lower = 0)
  if (!is.null(fisher) && !is.function(fisher)) {
    gdpar_abort("'fisher' must be a function of theta or NULL.",
                class = "gdpar_input_error")
  }
  th <- gdpar_geometry_thresholds()
  ladder <- .gdpar_orch_ladder()
  inform <- function(msg) if (isTRUE(verbose))
    gdpar_inform(msg, class = "gdpar_optin_message")

  # Build the engine sampling target from `target` when not supplied explicitly.
  if (is.null(geom_target)) {
    geom_target <- .gdpar_orch_derive_geom_target(target)
  } else if (!inherits(geom_target, "gdpar_geom_target")) {
    geom_target <- gdpar_geom_target(geom_target)
  }
  d <- geom_target$dim
  reference <- if (is.null(reference)) rep(0, d) else as.numeric(reference)
  if (length(reference) != d) {
    gdpar_abort(sprintf("'reference' must have length %d.", d),
                class = "gdpar_input_error")
  }

  inform("gdpar_geom_orchestrate is an opt-in geometry-adaptive controller; it does not modify any fit.")

  # --- initial diagnosis ---------------------------------------------------
  diag0 <- gdpar_geometry_diagnostic(target, n_grid = n_grid,
                                     seed = .gdpar_orch_seed(seed, 0L, 0L),
                                     verbose = FALSE, ...)
  sel <- .gdpar_orch_select_entry(diag0, level_map, entry_level, th)

  reproducibility <- list(
    seed = seed, gdpar_version = as.character(utils::packageVersion("gdpar")),
    budget = budget, criteria = criteria,
    initial_pathology = diag0$pathology, dim = d)

  # Out-of-scope diagnosis: short-circuit to a certificate naming the remedy.
  if (is.na(sel$key)) {
    inform(sprintf("Diagnosed '%s' is outside the geometry sampler ladder; emitting a certificate naming the proper remedy.",
                   diag0$pathology))
    presc <- .gdpar_orch_prescription(diag0, list(), is.function(fisher),
                                      budget, -Inf, character(0))
    cert <- .gdpar_orch_certificate("out_of_scope", diag0, list(),
                                    diag0$culprit, presc, reproducibility,
                                    list(fits = 0L, seconds = 0))
    obj <- list(status = "out_of_scope", certificate = cert, ledger = list(),
                diagnosis = diag0, best = NULL,
                budget_spent = list(fits = 0L, seconds = 0),
                reproducibility = reproducibility)
    class(obj) <- c("gdpar_geom_orchestration", "list")
    if (isTRUE(laplace_fallback)) {
      obj <- .gdpar_orch_attach_laplace(obj, geom_target, reference,
                                        laplace_draws, seed, inform)
    }
    return(obj)
  }

  # --- the armoured closed loop --------------------------------------------
  ladder_idx <- function(key) ladder[[key]]$index
  selected <- sel$key
  floor_index <- -Inf            # monotone ratchet.
  visited <- character(0)        # the (level, tier) memo.
  ledger <- list()
  fit_seconds <- numeric(0)
  fits <- 0L
  stall <- 0L
  best <- list(score = -Inf, draws = NULL, fit = NULL, key = NA_character_)
  warmup_draws <- NULL
  current_diag <- diag0

  run_fit <- function(metric, warmup, iter, slot, round, eps) {
    cap <- budget$max_seconds_per_fit
    if (is.finite(cap)) {
      setTimeLimit(elapsed = cap, transient = TRUE)
      on.exit(setTimeLimit(), add = TRUE)
    }
    t0 <- Sys.time()
    fit <- tryCatch(
      gdpar_geom_hmc(geom_target, metric = metric, epsilon = eps, L = budget$L,
                     n_iter = iter, n_warmup = warmup,
                     seed = .gdpar_orch_seed(seed, round, slot)),
      error = function(e)
        structure(list(error = conditionMessage(e)),
                  class = "gdpar_orch_fit_error"))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(fit = fit, elapsed = elapsed,
         failed = inherits(fit, "gdpar_orch_fit_error"))
  }

  for (round in seq_len(budget$max_rounds)) {
    # Budget admission: refuse a round we cannot afford (cost-aware).
    predicted <- if (length(fit_seconds)) mean(fit_seconds) else 0
    seconds_spent <- sum(fit_seconds)
    if (fits >= budget$max_fits ||
        (is.finite(budget$max_seconds) &&
         seconds_spent + predicted > budget$max_seconds) ||
        length(unique(vapply(ledger, function(r) r$level, character(1)))) >=
          budget$max_levels && !selected %in%
          vapply(ledger, function(r) r$level, character(1))) {
      inform("Budget admission stop: insufficient remaining budget for another level.")
      break
    }
    tier_key <- paste(selected, "probe", sep = ":")
    if (tier_key %in% visited) {
      # Cycle guard: climb the ladder instead of repeating a state.
      nxt <- .gdpar_orch_next_above(max(floor_index, ladder_idx(selected)))
      if (is.na(nxt)) { inform("Ladder exhausted (cycle guard)."); break }
      selected <- nxt
      next
    }
    visited <- c(visited, tier_key)

    # Per-level guard: build the metric; a malformed metric is recorded, not a
    # crash, and the controller escalates.
    mb <- .gdpar_orch_build_metric(selected, geom_target, fisher, reference,
                                   speed, rest_mass, warmup_draws)
    if (!isTRUE(mb$ok)) {
      inform(sprintf("Level '%s' unavailable: %s. Escalating.",
                     selected, mb$reason))
      ledger[[length(ledger) + 1L]] <- list(
        round = round, level = selected, level_index = ladder_idx(selected),
        diagnosis = current_diag$pathology, metric_ok = FALSE,
        reason = mb$reason, gate = NULL, accept = NA_real_,
        divergent_rate = NA_real_, ebfmi = NA_real_, seconds = 0,
        conflict = if (round == 1L) sel$conflict else NA)
      floor_index <- max(floor_index, ladder_idx(selected))
      nxt <- .gdpar_orch_next_above(floor_index)
      if (is.na(nxt)) break
      selected <- nxt
      next
    }

    # Coarse step-size search so the level is not failed for a mis-set step.
    eps_use <- if (!is.null(mb$metric$suggested_epsilon))
      mb$metric$suggested_epsilon else budget$epsilon
    if (isTRUE(budget$tune_epsilon)) {
      tu <- .gdpar_orch_tune_epsilon(geom_target, mb$metric, budget, seed, round)
      eps_use <- tu$epsilon
      fits <- fits + tu$fits; fit_seconds <- c(fit_seconds, tu$seconds)
    }
    # Cheap probe first (successive-halving spirit).
    pr <- run_fit(mb$metric, budget$probe_warmup, budget$probe_iter, 1L, round,
                  eps_use)
    fits <- fits + 1L; fit_seconds <- c(fit_seconds, pr$elapsed)
    probe_total <- budget$probe_warmup + budget$probe_iter
    probe_gate <- if (pr$failed) list(pass = FALSE, reasons = pr$fit$error)
                  else .gdpar_orch_success_gate(pr$fit, probe_total, criteria, 0)
    if (!pr$failed) {
      sc <- .gdpar_orch_score(pr$fit, probe_total)
      if (sc > best$score) {
        last <- pr$fit$draws[nrow(pr$fit$draws), ]
        best <- list(score = sc, draws = pr$fit$draws, fit = pr$fit,
                     key = selected)
        if (all(is.finite(last))) { reference <- last; warmup_draws <- pr$fit$draws }
      }
    }

    full_fit <- NULL; full_gate <- NULL; full_total <- NA_integer_
    if (isTRUE(probe_gate$pass)) {
      # The level earns the full budget; the final phase uses the frozen level.
      fl <- run_fit(mb$metric, budget$full_warmup, budget$full_iter, 2L, round,
                    eps_use)
      fits <- fits + 1L; fit_seconds <- c(fit_seconds, fl$elapsed)
      full_total <- budget$full_warmup + budget$full_iter
      full_gate <- if (fl$failed) list(pass = FALSE, reasons = fl$fit$error)
                   else .gdpar_orch_success_gate(fl$fit, full_total, criteria,
                                                 budget$hysteresis)
      full_fit <- fl$fit
      if (!fl$failed) {
        sc <- .gdpar_orch_score(fl$fit, full_total)
        if (sc > best$score) best <- list(score = sc, draws = fl$fit$draws,
                                          fit = fl$fit, key = selected)
      }
    }

    ledger[[length(ledger) + 1L]] <- list(
      round = round, level = selected, level_index = ladder_idx(selected),
      diagnosis = current_diag$pathology, metric_ok = TRUE,
      probe_gate = probe_gate, full_gate = full_gate,
      accept = if (!is.null(full_fit)) full_fit$accept_rate else
        if (!pr$failed) pr$fit$accept_rate else NA_real_,
      ebfmi = if (!is.null(full_fit)) full_fit$ebfmi else
        if (!pr$failed) pr$fit$ebfmi else NA_real_,
      seconds = pr$elapsed + (if (!is.null(full_fit)) fit_seconds[length(fit_seconds)] else 0),
      conflict = if (round == 1L) sel$conflict else NA)

    # Atomic checkpoint of the ledger for durability.
    if (!is.null(checkpoint_dir)) {
      .gdpar_orch_checkpoint(checkpoint_dir, ledger)
    }

    if (isTRUE(full_gate$pass)) {
      inform(sprintf("Resolved at level '%s'.", selected))
      obj <- list(
        status = "resolved", level = selected,
        level_index = ladder_idx(selected), metric = mb$metric,
        draws = full_fit$draws, fit = full_fit, gate = full_gate,
        ledger = ledger, diagnosis = diag0, best = best,
        budget_spent = list(fits = fits, seconds = sum(fit_seconds)),
        reproducibility = reproducibility)
      class(obj) <- c("gdpar_geom_orchestration", "list")
      return(obj)
    }

    # Failure: ratchet, then re-diagnose to choose the next level.
    floor_index <- max(floor_index, ladder_idx(selected))
    if (ladder_idx(selected) >=
        max(vapply(ladder[.gdpar_orch_ladder_keys()], function(l) l$index,
                   integer(1)))) {
      inform("Top of the ladder reached without success.")
      break
    }
    # Re-diagnose with a fresh pilot (warm-started reference already updated);
    # n_rediagnose >= 2 enables the stability check.
    rediag <- tryCatch(
      gdpar_geometry_diagnostic(target, n_grid = n_grid,
                                seed = .gdpar_orch_seed(seed, round, 9L),
                                verbose = FALSE, ...),
      error = function(e) NULL)
    safe_step <- .gdpar_orch_next_above(floor_index)
    if (is.null(rediag)) {
      next_sel <- safe_step
    } else {
      current_diag <- rediag
      stable <- TRUE
      if (budget$n_rediagnose >= 2L) {
        rediag2 <- tryCatch(
          gdpar_geometry_diagnostic(target, n_grid = n_grid,
                                    seed = .gdpar_orch_seed(seed, round, 10L),
                                    verbose = FALSE, ...),
          error = function(e) NULL)
        stable <- !is.null(rediag2) &&
          identical(rediag2$pathology, rediag$pathology)
      }
      if (!stable) {
        next_sel <- safe_step          # unstable re-diagnosis -> safe ladder step.
      } else {
        cand <- .gdpar_orch_pathology_to_key(rediag$pathology)
        if (is.na(cand) || ladder_idx(cand) <= floor_index) {
          next_sel <- safe_step        # ratchet: never below the floor.
        } else {
          next_sel <- cand             # the re-diagnosis-guided (possibly lateral) move.
        }
      }
    }
    if (is.na(next_sel)) { inform("No higher level to escalate to."); break }
    if (identical(next_sel, selected)) {
      stall <- stall + 1L
      if (stall >= budget$stall_limit) {
        inform("No-progress stall limit reached."); break
      }
      next_sel <- safe_step
      if (is.na(next_sel)) break
    } else {
      stall <- 0L
    }
    selected <- next_sel
  }

  # --- the budget is spent without success: certify the limit --------------
  tried_keys <- unique(vapply(ledger, function(r) r$level, character(1)))
  presc <- .gdpar_orch_prescription(current_diag, ledger, is.function(fisher),
                                    budget, best$score, tried_keys)
  cert <- .gdpar_orch_certificate("intractable_within_budget", current_diag,
                                  ledger, current_diag$culprit, presc,
                                  reproducibility,
                                  list(fits = fits, seconds = sum(fit_seconds)))
  inform("Certified limit: no level met the success gate within the budget.")
  obj <- list(status = "certified_limit", certificate = cert, ledger = ledger,
              diagnosis = diag0, best = best,
              budget_spent = list(fits = fits, seconds = sum(fit_seconds)),
              reproducibility = reproducibility)
  class(obj) <- c("gdpar_geom_orchestration", "list")
  if (isTRUE(laplace_fallback)) {
    obj <- .gdpar_orch_attach_laplace(obj, geom_target, reference,
                                      laplace_draws, seed, inform)
  }
  obj
}

# Derive the engine sampling target from the diagnostic target.
.gdpar_orch_derive_geom_target <- function(target) {
  if (inherits(target, "gdpar_geometry_target")) {
    inst <- target$make(target$default_n, target$default_difficulty)
    return(gdpar_geom_target(inst))
  }
  if (is.list(target) && (inherits(target$model, "CmdStanModel") ||
                          !is.null(target$model))) {
    if (is.null(target$dim)) {
      gdpar_abort(
        "A cmdstan model target needs 'dim' (the unconstrained dimension) to derive the engine target.",
        class = "gdpar_input_error")
    }
    return(gdpar_geom_target(object = target$model, dim = target$dim,
                             data = target$data))
  }
  gdpar_abort(
    paste("Cannot derive an engine target from 'target'; supply 'geom_target'",
          "explicitly (a gdpar_geom_target or a closure)."),
    class = "gdpar_input_error")
}

# Atomic ledger checkpoint (temp file then rename).
.gdpar_orch_checkpoint <- function(dir, ledger) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dir, fileext = ".rds")
  saveRDS(ledger, tmp)
  file.rename(tmp, file.path(dir, "orchestrate_ledger.rds"))
  invisible(NULL)
}

#' Print method for gdpar_geom_orchestration objects
#'
#' @param x A \code{gdpar_geom_orchestration}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_orchestration <- function(x, ...) {
  cat("<gdpar_geom_orchestration> status: ", x$status, "\n", sep = "")
  if (identical(x$status, "resolved")) {
    cat("  resolved at level: ", x$level, " (index ", x$level_index, ")\n",
        sep = "")
    cat("  draws: ", nrow(x$draws), " x ", ncol(x$draws),
        " | accept ", format(x$fit$accept_rate, digits = 3),
        " | E-BFMI ", format(x$fit$ebfmi, digits = 3), "\n", sep = "")
  } else {
    cat("  verdict: ", x$certificate$verdict, "\n", sep = "")
    cat("  levels attempted: ",
        paste(unique(vapply(x$ledger, function(r) r$level, character(1))),
              collapse = " -> "), "\n", sep = "")
    cat("  prescription items: ", length(x$certificate$prescription), "\n",
        sep = "")
  }
  if (!is.null(x$laplace)) {
    cat("  laplace fallback: '", x$laplace$fit_quality_label,
        "' (cond ", format(x$laplace$cond_unfloored, digits = 3),
        ", Pareto-k ", format(x$laplace$fit_quality$pareto_k, digits = 3),
        ")\n", sep = "")
  }
  cat("  budget spent: ", x$budget_spent$fits, " fits | ",
      format(x$budget_spent$seconds, digits = 3), " s\n", sep = "")
  invisible(x)
}

#' The certified-limit object of the geometry-adaptive orchestrator
#'
#' The first-class certificate emitted by \code{\link{gdpar_geom_orchestrate}}
#' when the budget is exhausted without a level meeting the success gate, or when
#' the diagnosis points outside the geometry sampler ladder. It separates the
#' demonstrated evidence in three rigour layers (\code{algebraic} = the geometry,
#' \code{statistical} = the per-level sampler diagnostics, \code{numerical} = the
#' budget and fits) from a \code{prescription} of conjectured, falsifiable fixes,
#' and carries a \code{reproducibility} block. This print method summarises it.
#'
#' @param x A \code{gdpar_geom_certificate}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{gdpar_geom_orchestrate}}.
#' @aliases gdpar_geom_certificate
#' @export
print.gdpar_geom_certificate <- function(x, ...) {
  cat("<gdpar_geom_certificate> verdict: ", x$verdict, "\n", sep = "")
  cat("  [algebraic] pathology: ", x$algebraic$pathology,
      " | condition ", format(x$algebraic$condition_number, digits = 3),
      " | grows with n: ", x$algebraic$grows_with_n, "\n", sep = "")
  cat("  [statistical] levels recorded: ", length(x$statistical), "\n", sep = "")
  cat("  [numerical] fits: ", x$numerical$fits, " | seconds: ",
      format(x$numerical$seconds, digits = 3), "\n", sep = "")
  cat("  [prescription] (conjectured):\n")
  for (p in x$prescription) {
    cat("    - ", p$action, "\n", sep = "")
  }
  invisible(x)
}

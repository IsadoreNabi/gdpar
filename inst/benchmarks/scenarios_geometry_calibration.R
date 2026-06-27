## Block RG, RG.1.c: calibration scenario grid for the posterior-geometry
## diagnostic (gdpar_geometry_diagnostic + gdpar_geometry_suite).
##
## The expensive work of RG.1.c is COLLECTING the size-invariant signals by
## running cheap NUTS pilots through cmdstan; the threshold sweep itself is a
## cheap, pure function of the stored signals (see rg_calibration_analysis.R).
## This file defines the three calibration axes and the per-cell runner; it is
## a benchmark script (sourced under devtools::load_all), not package code, so
## the default fit path stays bit-identical.
##
## Three axes (Q2 = "mix 3+4 adaptive", session B9.23):
##   1. difficulty: a per-target severity ladder (low/mid/high pathology
##      intensity, mid ~ the suite default_difficulty). This drives the
##      error-vs-severity curve and the honesty caveat (a classifier may only
##      catch SEVERE pathologies -- e.g. a deep funnel but not a shallow one).
##   2. pilot budget: cheap (100/100) vs default (150/150) warmup/sampling, to
##      report how diagnostic accuracy degrades as the pilots get cheaper (the
##      cost-vs-accuracy curve of the diagnostic itself).
##   3. replicas: distinct seeds per cell. Allocated ADAPTIVELY (see the driver
##      and analysis): Stage 1 gives every cell R0 replicas; Stage 2 tops up the
##      cells whose Wilson confidence interval for the correct-classification
##      rate is WIDEST (the p ~ 0.5 decision-boundary band) -- a minimax-robust
##      allocation that spends the budget where the error rate is least certain.

## Per-target severity ladder. The severity RANK (low/mid/high), not the raw
## knob, indexes the error-vs-severity curve, because the knob direction varies
## by target: larger is more severe for the funnel (scale_v), anisotropy
## (kappa), flat direction (tau) and multimodal separation (sep); SMALLER nu is
## heavier-tailed; SMALLER alpha pins the boundary harder; G0/G4 sweep the
## dimension (G0 is the negative control, G4's signature is the n-growth curve,
## not the difficulty knob).
rg_difficulty_ladder <- function() {
  list(
    G0_isotropic           = c(low = 3,   mid = 5,   high = 10),
    G1_anisotropic         = c(low = 20,  mid = 100, high = 500),
    G2_funnel              = c(low = 1,   mid = 3,   high = 6),
    G3_heavy_tails         = c(low = 4,   mid = 2,   high = 1.3),
    G4_quasi_deterministic = c(low = 3,   mid = 5,   high = 8),
    G5_multimodal          = c(low = 2.5, mid = 4,   high = 6),
    G6_boundary            = c(low = 0.6, mid = 0.3, high = 0.1),
    G7_flat_direction      = c(low = 30,  mid = 100, high = 300)
  )
}

## Pilot-budget levels (warmup, sampling). Chains stay at 4 (needed for the
## multimodality signal).
rg_pilot_budgets <- function() {
  list(
    cheap   = c(warmup = 100L, sampling = 100L),
    default = c(warmup = 150L, sampling = 150L)
  )
}

## The full cell-base grid (target x difficulty x pilot), independent of reps.
rg_cell_grid <- function() {
  ladder  <- rg_difficulty_ladder()
  budgets <- rg_pilot_budgets()
  rows <- list()
  for (tid in names(ladder)) {
    diffs <- ladder[[tid]]
    for (dl in names(diffs)) {
      for (pl in names(budgets)) {
        rows[[length(rows) + 1L]] <- data.frame(
          target_id   = tid,
          diff_level  = dl,
          diff_value  = unname(diffs[[dl]]),
          pilot_level = pl,
          warmup      = unname(budgets[[pl]]["warmup"]),
          sampling    = unname(budgets[[pl]]["sampling"]),
          cell_base   = sprintf("%s__%s__%s", tid, dl, pl),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

## Reproducible per-replica seed: a deterministic function of the cell base and
## the replica index, so reruns and resumes reproduce the same pilots.
rg_cell_seed <- function(cell_base, rep) {
  base <- sum(as.integer(charToRaw(cell_base)))
  as.integer(20260602L + base * 101L + rep * 7L)
}

## Run one calibration cell: a single gdpar_geometry_diagnostic at the cell's
## difficulty and pilot budget, with the replica seed. Returns the compact,
## storable object from which any threshold set can be re-evaluated offline
## (the full per-n signals data frame + the n-curve slope + the culprit + the
## ground-truth label).
rg_run_cell <- function(suite, crow, rep, seed) {
  tg <- suite[[crow$target_id]]
  t0 <- Sys.time()
  diag <- gdpar_geometry_diagnostic(
    tg,
    difficulty     = crow$diff_value,
    pilot_warmup   = as.integer(crow$warmup),
    pilot_sampling = as.integer(crow$sampling),
    chains         = 4L,
    seed           = seed,
    verbose        = FALSE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(
    cell_id      = sprintf("%s__rep%02d", crow$cell_base, rep),
    cell_base    = crow$cell_base,
    target_id    = crow$target_id,
    diff_level   = crow$diff_level,
    diff_value   = crow$diff_value,
    pilot_level  = crow$pilot_level,
    warmup       = as.integer(crow$warmup),
    sampling     = as.integer(crow$sampling),
    rep          = rep,
    seed         = seed,
    signals      = diag$signals,
    slope        = diag$difficulty_curve$slope,
    culprit      = diag$culprit,
    ground_truth = diag$ground_truth$pathology,
    default_pathology = diag$pathology,
    default_correct   = isTRUE(diag$correct),
    elapsed      = elapsed,
    status       = if (all(diag$signals$failed)) "all_pilots_failed"
                   else if (any(diag$signals$failed)) "some_pilots_failed"
                   else "ok"
  )
}

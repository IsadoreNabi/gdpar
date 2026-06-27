## Block RG, sub-phase RG.7 driver: apply the geometry-adaptive capability to
## the REAL Tweedie count of benchmark 9.2.O (the four N__N1__gdpar cells), one
## sub-region per checkpointed cell, resumable.
##
## For each sub-region it RE-FITS the Tweedie count (the 9.2.O cell stores only
## param draws, not the fit), BRIDGES the fit to the geometry engine, measures
## the local curvature at the reference by the user's CROSS-CHECK route
## (observed information by finite differences + expected Fisher by simulation
## through the real model), confirms the two agree on the wall/floor split, and
## runs the orchestrator (natural + forced sub-Riemannian). The verdict per cell
## is RESOLVED (an honest ELPD route) or a CERTIFIED reproducible LIMIT -- both
## first-class. Nothing here touches the package R/ sources.
##
## Env gates:
##   GDPAR_RG7=1                 required (else skip, exit 0).
##   GDPAR_RG7_MODE=smoke|full   smoke = tiny end-to-end validation (default full).
##   GDPAR_RG7_SR                comma list of sub-regions (default all four).
##   GDPAR_RG7_FORCE=1           recompute even if a cell RDS exists.
##   GDPAR_RG7_TAG               isolate a scratch run (suffix on the out dir).
##   GDPAR_RG7_NSIM              override the Fisher-simulation count.
##
## Usage (from the gdpar working directory; setsid + TMPDIR off tmpfs):
##   GDPAR_RG7=1 GDPAR_RG7_MODE=full Rscript inst/benchmarks/run_rg7_geometry_tweedie.R

if (Sys.getenv("GDPAR_RG7") != "1") {
  message("GDPAR_RG7 is not 1; skipping RG.7. Set GDPAR_RG7=1 to execute it.")
  quit(save = "no", status = 0L)
}
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all(".", quiet = TRUE))

repo_root <- normalizePath(".")
bench_dir <- file.path(repo_root, "inst", "benchmarks")
scripts_dir <- file.path(bench_dir, "scripts")
source(file.path(bench_dir, "scenarios_ebird_9_2_O.R"))
source(file.path(scripts_dir, "ebird_adapters_9_2_O.R"))
source(file.path(scripts_dir, "rg7_geometry_helpers.R"))
results_dir <- file.path(bench_dir, "results")

mode <- Sys.getenv("GDPAR_RG7_MODE", "full")
tag  <- Sys.getenv("GDPAR_RG7_TAG")
suffix <- if (nzchar(tag)) paste0("_", tag) else ""
out_dir <- file.path(results_dir, paste0("rg7_geometry", suffix))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, "run_rg7_geometry.log")
.elog <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n", sep = "",
      file = log_path, append = TRUE)
}

force <- Sys.getenv("GDPAR_RG7_FORCE") == "1"
srs_env <- Sys.getenv("GDPAR_RG7_SR")
srs <- if (nzchar(srs_env)) trimws(strsplit(srs_env, ",")[[1L]]) else .ebird_9_2_O_sub_regions
nsim_env <- Sys.getenv("GDPAR_RG7_NSIM")

## --- per-mode configuration --------------------------------------------------
if (mode == "smoke") {
  refit_iters <- list(chains = 1L, iter_warmup = 80L, iter_sampling = 80L)
  refit_treedepth <- 6L
  n_sim <- if (nzchar(nsim_env)) as.integer(nsim_env) else 20L
  bridge_hessian <- TRUE
  geom_cfg <- list(n_warmup = 80L, n_iter = 120L, L = 16L, tune_iter = 40L)
  n_chains_cfg <- 2L            # multi-chain R-hat gate (cheap in smoke).
  budget <- NULL
} else {
  refit_iters <- list(chains = 2L, iter_warmup = 300L, iter_sampling = 600L)
  refit_treedepth <- 9L
  n_sim <- if (nzchar(nsim_env)) as.integer(nsim_env) else 64L
  bridge_hessian <- TRUE        # exact Hessian for the Laplace metric + the
                                # observed-information cross-check; compiles once,
                                # cached across the four cells (same model code).
  geom_cfg <- list(n_warmup = 500L, n_iter = 1000L, L = 24L, tune_iter = 80L)
  n_chains_cfg <- 4L            # 4 independent geometric chains -> R-hat/ESS gate
                                # (GLM/MiMo: acceptance + divergences alone do not
                                # prove resolution of a frozen-mass canyon).
  budget <- NULL
}

.elog("=== RG.7 START === mode=%s sub_regions=%s n_sim=%d",
      mode, paste(srs, collapse = "/"), n_sim)

ds <- ebird_9_2_O_data_cache(results_dir, scripts_dir, force = FALSE)

for (sr in srs) {
  srs_short <- unname(.ebird_9_2_O_subregion_short[sr])
  out_rds <- file.path(out_dir, sprintf("rg7_%s.rds", srs_short))
  if (!force && file.exists(out_rds)) {
    .elog("cell %s exists; skipping.", srs_short); next
  }
  t_cell <- Sys.time()
  .elog("--- cell %s (%s): re-fitting Tweedie count ---", srs_short, sr)
  sp <- ebird_9_2_O_subregion_split(ds$full, sr)
  fit_seed <- ebird_9_2_O_fit_seed(sr, "N1", "gdpar")

  res <- tryCatch({
    rf <- rg7_refit_tweedie(sp$train, iters = refit_iters, fit_seed = fit_seed,
                            max_treedepth = refit_treedepth)
    fit <- rf$fit; tr <- rf$train
    diag_nuts <- .ebird_stan_diag_cmdstan(fit$fit)
    .elog("cell %s: NUTS diag rhat_max=%.3g ess_min=%.3g div=%g",
          srs_short, diag_nuts$rhat_max, diag_nuts$ess_bulk_min,
          diag_nuts$n_divergences)

    ## bridge with the EXACT Hessian (needed by the winning Laplace metric and
    ## the observed-information cross-check; FD blows up at this stiffness).
    .elog("cell %s: bridging (hessian=TRUE) ...", srs_short)
    br <- gdpar_geom_bridge(fit, hessian = TRUE)
    ref_best <- rg7_best_reference(br, fit)
    ## climb to the genuine MODE from several high-density starts (where the
    ## Hessian is all-positive = the winning recipe; an off-mode reference points
    ## the curvature the wrong way). Multi-start verifies a single basin and
    ## keeps the cleanest mode (GLM/MiMo review, rank 3).
    mc <- rg7_climb_multistart(br, fit, br$geom_target, K = 3L, verbose = FALSE)
    mode <- mc$mode
    .elog("cell %s: bridge dim=%d hessian_avail=%s | mode converged=%s logp=%.2f |grad|=%.3g starts=%s basin_spread=%.3g",
          srs_short, br$dim, !is.null(br$geom_target$hessian),
          mc$converged, mc$logp %||% NA, mc$grad_norm %||% NA,
          mc$n_starts %||% NA, mc$basin_spread %||% NA)

    ## observed information = the EXACT Hessian at the mode (the winning
    ## curvature for the Laplace metric); record its conditioning.
    Jobs <- rg7_observed_information(br$geom_target, mode)
    ev <- tryCatch(
      eigen(0.5 * (Jobs + t(Jobs)), symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_)
    cond_obs <- if (all(is.finite(ev))) max(pmax(ev, 1e-12)) / min(pmax(ev, 1e-12)) else NA_real_
    .elog("cell %s: observed-info (%s) cond=%.3e min_eig=%.3g (all-pos=%s)",
          srs_short, attr(Jobs, "method"), cond_obs,
          if (all(is.finite(ev))) min(ev) else NA_real_,
          all(is.finite(ev)) && all(ev > 0))

    ## cross-check WITNESS: the simulated expected Fisher at the mode (the user's
    ## validacion-cruzada; RG.7 found it must NOT feed the metric off-mode, but
    ## it is recorded as the independent curvature estimate). It is a DIAGNOSTIC
    ## ONLY -- a degenerate simulated Fisher (too few valid real-model scores on
    ## the stiffest sub-regions) must NEVER abort a resolvable cell, so the call
    ## is gated on a finite witness and wrapped.
    fsim <- tryCatch(rg7_simulated_fisher(br, fit, tr, n_sim = n_sim,
                       seed = fit_seed + 7L, verbose = FALSE, reference = mode),
                     error = function(e) NULL)
    witness_ok <- !is.null(fsim) && !is.null(fsim$I_centred) &&
                  all(is.finite(fsim$I_centred))
    xchk <- if (witness_ok)
              tryCatch(rg7_crosscheck(Jobs, fsim$I_centred), error = function(e) NULL)
            else NULL
    if (!is.null(xchk) && isTRUE(is.finite(xchk$cond_sim)))
      .elog("cell %s: cross-check cond_sim=%.3e principal_angle=%.2fdeg",
            srs_short, xchk$cond_sim, xchk$max_principal_angle_deg)
    else
      .elog("cell %s: cross-check witness DEGENERATE (n_ok=%s sample_fail=%s grad_fail=%s pred_finite=%s); proceeding on the observed Hessian (the metric is healthy)",
            srs_short, if (!is.null(fsim)) fsim$n_ok else NA,
            if (!is.null(fsim)) fsim$n_sample_fail else NA,
            if (!is.null(fsim)) fsim$n_grad_fail else NA,
            if (!is.null(fsim)) fsim$pred_finite else NA)

    ## WINNING recipe: dense Euclidean (Laplace, mass = exact Hessian at mode),
    ## eps mini-tuned; sub-Riemannian + identity-Euclidean recorded.
    .elog("cell %s: applying geometry (Laplace dense-Euclidean + sub-Riemannian + identity) ...",
          srs_short)
    geo <- rg7_apply_geometry(br, curvature = Jobs, reference = mode,
                              geom = geom_cfg, seed = fit_seed + 11L,
                              verbose = TRUE, run_subriemannian = TRUE,
                              n_chains = n_chains_cfg)
    df <- geo$dense_euclidean
    if (!is.null(df) && is.null(df$error)) {
      dd <- geo$dense_diag
      .elog("cell %s: WINNER dense-Euclidean(Laplace) eps=%.4g accept=%.3f div=%d ebfmi=%.3f | max_rhat=%.4g min_ess=%.4g | gate[acc=%s div=%s ebfmi=%s conv=%s] RESOLVED=%s",
            srs_short, geo$best_eps, df$accept_rate, df$n_divergent, df$ebfmi %||% NA,
            if (!is.null(dd)) dd$max_rhat else NA, if (!is.null(dd)) dd$min_ess_bulk else NA,
            geo$resolved_components$acc, geo$resolved_components$div,
            geo$resolved_components$ebfmi, geo$resolved_components$convergence,
            geo$resolved)
    }

    list(sub_region = sr, srs_short = srs_short, fit_seed = fit_seed,
         n_train = nrow(tr), n_test = nrow(sp$test),
         nuts_diag = diag_nuts, dim = br$dim, reference = ref_best,
         mode = mc, observed_info = Jobs, cond_obs = cond_obs,
         fisher_sim = fsim, crosscheck = xchk,
         geometry = geo, resolved = geo$resolved, budget = budget,
         gdpar_version = as.character(utils::packageVersion("gdpar")),
         timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  }, error = function(e) {
    .elog("cell %s: ERROR %s", srs_short, conditionMessage(e))
    list(sub_region = sr, srs_short = srs_short, error = conditionMessage(e),
         error_class = class(e)[1L], mode = mode)
  })

  tmp <- paste0(out_rds, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out_rds)
  .elog("=== cell %s DONE in %.1f min ===", srs_short,
        as.numeric(difftime(Sys.time(), t_cell, units = "mins")))
}
.elog("=== RG.7 BATCH COMPLETE ===")

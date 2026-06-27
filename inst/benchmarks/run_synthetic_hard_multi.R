## Block 6 runner: synthetic adversarial validation TOP4 with
## coord-wise factorization.
##
## Iterates over the 9 scenarios x 5 competitors (gdpar, brms, mgcv,
## INLA, rstanarm), with brms running a second pass on S8 with
## set_rescor(TRUE) under the method label "brms_rescor". Each
## (scenario, competitor) cell yields a per-coordinate metric block
## (RMSE, MAE, CRPS for gaussian; AUC, log_loss, Brier for
## bernoulli), one aggregated elpd_loo row with Pareto-k diagnostics,
## one wall-time row, and competitor-specific cross-check rows
## (AIC for mgcv, elpd_cpo for INLA).
##
## Outputs: a single long-tidy CSV at
##   inst/benchmarks/results/synthetic_hard_multi_results.csv
##
## Env gate: GDPAR_BENCH_MULTI=1
## Wall budget: ~1-2 hours total on a Ryzen 9 5980HX.
##
## Usage (from the gdpar working directory):
##   GDPAR_BENCH_MULTI=1 Rscript inst/benchmarks/run_synthetic_hard_multi.R
##
## To run a subset of scenarios:
##   GDPAR_BENCH_MULTI=1 \
##     GDPAR_BENCH_SCENARIOS="S0,S1" \
##     Rscript inst/benchmarks/run_synthetic_hard_multi.R

if (Sys.getenv("GDPAR_BENCH_MULTI") != "1") {
  message(
    "GDPAR_BENCH_MULTI is not set to 1; skipping the synthetic ",
    "multi bench run. Set GDPAR_BENCH_MULTI=1 to execute it."
  )
  quit(save = "no", status = 0L)
}

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to load gdpar via load_all().")
}
suppressMessages(devtools::load_all(".", quiet = TRUE))

repo_root <- normalizePath(".")
scripts_dir <- file.path(repo_root, "inst", "benchmarks", "scripts")
source(file.path(scripts_dir, "bench_multi_helpers.R"))
source(file.path(scripts_dir, "bench_competitor_gdpar.R"))
source(file.path(scripts_dir, "bench_competitor_brms.R"))
source(file.path(scripts_dir, "bench_competitor_mgcv.R"))
source(file.path(scripts_dir, "bench_competitor_inla.R"))
source(file.path(scripts_dir, "bench_competitor_rstanarm.R"))
source(file.path(repo_root, "inst", "benchmarks",
                  "synthetic_hard_multi.R"))

results_dir <- file.path(repo_root, "inst", "benchmarks", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(results_dir, "synthetic_hard_multi_results.csv")

subset_env <- Sys.getenv("GDPAR_BENCH_SCENARIOS")
all_gens <- scenarios_synthetic_hard_multi()
if (nzchar(subset_env)) {
  wanted <- strsplit(subset_env, ",", fixed = TRUE)[[1L]]
  wanted <- trimws(wanted)
  all_gens <- all_gens[wanted]
  if (length(all_gens) == 0L) {
    stop(sprintf(
      "GDPAR_BENCH_SCENARIOS = '%s' selected no scenarios; valid ids: %s",
      subset_env,
      paste(names(scenarios_synthetic_hard_multi()), collapse = ", ")
    ))
  }
}

competitors <- list(
  gdpar    = list(fit = bench_gdpar_fit,
                   metrics = bench_gdpar_metrics),
  brms     = list(fit = bench_brms_fit,
                   metrics = bench_brms_metrics),
  mgcv     = list(fit = bench_mgcv_fit,
                   metrics = bench_mgcv_metrics),
  inla     = list(fit = bench_inla_fit,
                   metrics = bench_inla_metrics),
  rstanarm = list(fit = bench_rstanarm_fit,
                   metrics = bench_rstanarm_metrics)
)

run_one_competitor <- function(name, fit_fn, metrics_fn, data,
                                family_spec, scenario_meta,
                                extra_args = list()) {
  message(sprintf("  [%s] fitting ...", name))
  res <- tryCatch(
    do.call(fit_fn,
            c(list(data = data,
                    family_spec = family_spec,
                    scenario_meta = scenario_meta),
              extra_args)),
    error = function(e) {
      warning(sprintf("[%s] fit failed: %s", name, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) {
    return(bench_record(
      scenario_meta$scenario_id, scenario_meta$p, name,
      family_spec$brms_family, "fit_failed", "agg", NA_real_,
      convergence_flag = "fit_failed",
      extra = "competitor fit_fn raised an error"
    ))
  }
  message(sprintf("  [%s] computing metrics ...", name))
  rows <- tryCatch(
    metrics_fn(res, scenario_meta$.test_data, scenario_meta),
    error = function(e) {
      warning(sprintf("[%s] metrics failed: %s", name,
                       conditionMessage(e)))
      bench_record(
        scenario_meta$scenario_id, scenario_meta$p, name,
        family_spec$brms_family, "metrics_failed", "agg", NA_real_,
        convergence_flag = "metrics_failed",
        wall_seconds = res$wall_seconds,
        extra = conditionMessage(e))
    }
  )
  rows
}

t_start <- Sys.time()
all_rows <- list()
for (sc_name in names(all_gens)) {
  message(sprintf("\n=== Scenario %s ===", sc_name))
  gen <- all_gens[[sc_name]]
  sc <- gen()
  sc$scenario_meta$.test_data <- sc$test
  for (comp_name in names(competitors)) {
    comp <- competitors[[comp_name]]
    rows <- run_one_competitor(comp_name, comp$fit, comp$metrics,
                                 sc$train, sc$family_spec,
                                 sc$scenario_meta)
    if (!is.null(rows)) all_rows[[length(all_rows) + 1L]] <- rows
  }
  if (isTRUE(sc$scenario_meta$brms_rescor_recovery)) {
    message("  [brms_rescor] additional fit with set_rescor(TRUE) ...")
    rows <- run_one_competitor(
      "brms_rescor",
      function(data, family_spec, scenario_meta, ...) {
        bench_brms_fit(data, family_spec, scenario_meta, rescor = TRUE)
      },
      function(fit_obj, test, scenario_meta) {
        bench_brms_metrics(fit_obj, test, scenario_meta,
                            method_label = "brms_rescor")
      },
      sc$train, sc$family_spec, sc$scenario_meta)
    if (!is.null(rows)) all_rows[[length(all_rows) + 1L]] <- rows
  }
  rm(sc); gc(verbose = FALSE)
}
t_end <- Sys.time()

results <- do.call(rbind, all_rows)
write.csv(results, file = out_path, row.names = FALSE)

message(sprintf(
  "\nWrote %d rows to %s\nWall total: %.1f min",
  nrow(results), out_path,
  as.numeric(difftime(t_end, t_start, units = "mins"))
))

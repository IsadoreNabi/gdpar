## Block 7 runner: external validation eBird Status & Trends p = 2
## coord-wise factorization (logit-encounter, log-count) with state
## grouping (4 sub-regions NE-USA breeding x 5 forest passerines).
##
## Iterates over the single ebird_path_a_multi scenario x 5
## competitors (gdpar, brms, mgcv, INLA, rstanarm). Each
## (scenario, competitor) cell yields per-coordinate RMSE / MAE /
## CRPS rows plus an aggregated elpd_loo row with Pareto-k
## diagnostics, a wall-time row, and competitor-specific
## cross-checks.
##
## Outputs: a single long-tidy CSV at
##   inst/benchmarks/results/ebird_path_a_multi_results.csv
##
## Env gate: GDPAR_EBIRD_BENCH=1
## Wall budget: ~1-3 hours total on a Ryzen 9 5980HX.
##
## Usage (from the gdpar working directory):
##   GDPAR_EBIRD_BENCH=1 Rscript inst/benchmarks/run_ebird_a_multi.R

if (Sys.getenv("GDPAR_EBIRD_BENCH") != "1") {
  message(
    "GDPAR_EBIRD_BENCH is not set to 1; skipping the eBird external ",
    "validation run. Set GDPAR_EBIRD_BENCH=1 to execute it."
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
source(file.path(repo_root, "inst", "benchmarks", "synthetic_hard_multi.R"))
source(file.path(scripts_dir, "ebird_loader.R"))

results_dir <- file.path(repo_root, "inst", "benchmarks", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(results_dir, "ebird_path_a_multi_results.csv")

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

message("=== Building eBird path-A multi dataset ===")
ds <- build_ebird_path_a_multi_dataset()
ds$scenario_meta$.test_data <- ds$test

message(sprintf(
  "Train: %d rows, test: %d rows, groups: %d.",
  nrow(ds$train), nrow(ds$test), nlevels(ds$train$group)
))

t_start <- Sys.time()
all_rows <- list()
for (comp_name in names(competitors)) {
  comp <- competitors[[comp_name]]
  rows <- run_one_competitor(comp_name, comp$fit, comp$metrics,
                               ds$train, ds$family_spec,
                               ds$scenario_meta)
  if (!is.null(rows)) all_rows[[length(all_rows) + 1L]] <- rows
}
t_end <- Sys.time()

results <- do.call(rbind, all_rows)
write.csv(results, file = out_path, row.names = FALSE)

message(sprintf(
  "\nWrote %d rows to %s\nWall total: %.1f min",
  nrow(results), out_path,
  as.numeric(difftime(t_end, t_start, units = "mins"))
))

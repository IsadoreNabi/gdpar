# ============================================================
# gdpar / inst/benchmarks/full_validation.R
# Orchestrator that runs the full validation and benchmarking
# pipeline. Each section can be enabled or disabled via the
# RUN_FLAGS list at the top of the file.
#
# Recommended usage:
#   - Open this file in RStudio, edit RUN_FLAGS as needed
#     (the eBird section is OFF by default because it requires
#     the ebirdst data download), then select-all and press
#     Ctrl+Enter.
#   - Alternatively, set the working directory to either the
#     gdpar package root or to inst/benchmarks/, then
#     source("...").
#
# From the shell:
#   Rscript inst/benchmarks/full_validation.R
#
# Total runtime depends on which sections are enabled. Reasonable
# upper bounds on a modern multi-core CPU:
#   quick_smoke           5 min
#   synthetic_recovery   20 min
#   synthetic_hard       40 min
#   benchmark_competitors 60 min
#   ebird_validation     90 min  (plus initial download)
# ============================================================

RUN_FLAGS <- list(
  quick_smoke           = TRUE,
  synthetic_recovery    = TRUE,
  synthetic_hard        = TRUE,
  benchmark_competitors = TRUE,
  ebird_validation      = FALSE
)

bench_dir <- local({
  candidates <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                    error = function(e) "")
    if (nzchar(ctx) && file.exists(ctx)) {
      candidates <- c(candidates, dirname(ctx))
    }
  }
  sourced <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(sourced) && nzchar(sourced)) {
    candidates <- c(candidates,
                    dirname(normalizePath(sourced, mustWork = FALSE)))
  }
  candidates <- c(candidates,
                  getwd(),
                  file.path(getwd(), "inst", "benchmarks"))
  inst_dir <- system.file("benchmarks", package = "gdpar")
  if (nzchar(inst_dir)) candidates <- c(candidates, inst_dir)
  found <- ""
  for (d in candidates) {
    if (nzchar(d) && file.exists(file.path(d, "helpers.R"))) {
      found <- d; break
    }
  }
  if (!nzchar(found)) {
    stop("Could not locate the gdpar benchmarks directory. Set the ",
         "working directory to the package root or to inst/benchmarks/, ",
         "or install gdpar so that system.file() resolves it.")
  }
  found
})

source(file.path(bench_dir, "helpers.R"))

bench_section("gdpar full validation pipeline")
bench_status(sprintf("Run flags: %s",
                       paste(names(Filter(isTRUE, RUN_FLAGS)),
                             collapse = ", ")))
bench_status(sprintf("Benchmarks directory: %s", bench_dir))
bench_status(sprintf("Working directory: %s", getwd()))

if (!bench_require("gdpar")) stop("gdpar must be installed first.")

if (isTRUE(RUN_FLAGS$quick_smoke)) {
  source(file.path(bench_dir, "quick_smoke.R"), local = TRUE)
}

if (isTRUE(RUN_FLAGS$synthetic_recovery)) {
  source(file.path(bench_dir, "synthetic_recovery.R"), local = TRUE)
}

if (isTRUE(RUN_FLAGS$synthetic_hard)) {
  source(file.path(bench_dir, "synthetic_hard.R"), local = TRUE)
}

if (isTRUE(RUN_FLAGS$benchmark_competitors)) {
  source(file.path(bench_dir, "benchmark_competitors.R"), local = TRUE)
}

if (isTRUE(RUN_FLAGS$ebird_validation)) {
  source(file.path(bench_dir, "ebird_validation.R"), local = TRUE)
}

bench_section("Full validation pipeline complete")
cat("Results are stored in the global environment under:\n")
cat("  synthetic_recovery_results, synthetic_hard_results,\n")
cat("  competitor_results, ebird_results\n")
cat("(only the ones whose section was enabled).\n")

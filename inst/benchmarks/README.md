# gdpar benchmark and validation suite

This directory contains the standalone validation and benchmarking
scripts for the **gdpar** package. They are not part of the package
API; they exercise the API on synthetic and real data, compare
performance against canonical R competitors, and report results in a
human-readable form.

## Files

| File | Purpose | Approximate runtime |
|---|---|---|
| `helpers.R` | Shared utilities: section reporting, conditional package loading, synthetic data generators with controlled difficulty, predictive metrics (RMSE, MAE, CRPS, AUC, log-loss, Brier, PIT, calibration, ELPD-LOO via PSIS), coverage. | (sourced by all others) |
| `quick_smoke.R` | Fast end-to-end smoke test of every public function. | 2 to 5 minutes |
| `synthetic_recovery.R` | Parameter recovery on synthetic data with known truth. | 10 to 20 minutes |
| `synthetic_hard.R` | Adversarial synthetic data tests (heterogeneity, near-collinearity, heavy tails, zero-inflation). | 20 to 40 minutes |
| `benchmark_competitors.R` | Comparison against mgcv, rstanarm, brms (structural) and ranger, xgboost (predictive agnostic). | 30 to 60 minutes |
| `ebird_validation.R` | External validation on eBird Status & Trends abundance and occurrence for American Robin in NE USA. | 30 to 90 minutes (plus initial download) |
| `full_validation.R` | Orchestrator that runs the above in order, with each section toggleable via flags. | sum of above |

## Prerequisites

### Required

- R version 4.2.0 or above.
- The **gdpar** package installed (`devtools::install_local()` from the
  package root, or `R CMD INSTALL .`).
- The R packages declared in `Imports` of the `DESCRIPTION` file.

### Required for end-to-end fitting

- **cmdstanr** package and a working CmdStan installation:

```r
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           "https://cloud.r-project.org"))
cmdstanr::install_cmdstan(cores = 4)
```

### Required for the competitor benchmark (optional but recommended)

- `mgcv` (base R install in most distributions),
- `rstanarm`,
- `brms` (heavier; can be skipped),
- `ranger`,
- `xgboost`,
- `scoringRules` (for the CRPS metric).

```r
install.packages(c("rstanarm", "ranger", "xgboost", "scoringRules"),
                 repos = "https://cloud.r-project.org")
```

`brms` is heavier to install:

```r
install.packages("brms", repos = "https://cloud.r-project.org")
```

If a competitor package is not installed, the script reports a
`[skip]` line for that competitor and continues with the others.

### Required for the eBird validation (optional)

- The `ebirdst` and `terra` R packages.
- A valid eBird Status & Trends API key in the environment variable
  `EBIRDST_KEY`.
- Approximately 1 to 3 GB of disk space for the downloaded
  low-resolution data of one species.

```r
install.packages(c("ebirdst", "terra"),
                 repos = "https://cloud.r-project.org")
Sys.setenv(EBIRDST_KEY = "your_key_here")
```

## How to run

From the package root directory, in an R session:

```r
setwd("/path/to/gdpar")

source("inst/benchmarks/quick_smoke.R")
```

Or run the full pipeline in one call:

```r
source("inst/benchmarks/full_validation.R")
```

To enable or disable individual sections, edit the `RUN_FLAGS` list at
the top of `full_validation.R`. The `ebird_validation` flag is `FALSE`
by default because it requires the data download.

From the shell:

```bash
cd /path/to/gdpar
Rscript inst/benchmarks/full_validation.R
```

## Reading the output

Each script writes section banners and a final summary table to
standard output. End-of-section results are also assigned to the
global environment under named objects:

- `synthetic_recovery_results`,
- `synthetic_hard_results`,
- `competitor_results`,
- `ebird_results`.

To redirect the output to a log file:

```bash
Rscript inst/benchmarks/full_validation.R > validation_log.txt 2>&1
```

## Reproducibility notes

- Each generator in `helpers.R` accepts a `seed` argument and uses
  `set.seed(seed)` internally. The seeds are fixed by default so that
  multiple runs produce identical results.
- The recovery script (`synthetic_recovery.R`) uses `REPLICATES = 5`
  by default. Increase this number at the top of the script for
  tighter Monte Carlo bounds on the bias and RMSE estimates of the
  recovery experiments.
- The eBird script downloads data via `ebirdst::ebirdst_download_status`,
  which is idempotent across calls; subsequent runs reuse the local
  cache without redownloading.

## Methodological caveats

- All gdpar fits in this version use Path 1 (hierarchical Bayesian via
  cmdstanr) with finite-dimensional parametric AMM specifications and
  a single individual-specific parameter. The competitor benchmark
  comparisons reflect this scope.
- The competitor benchmark reports results in two clearly separated
  groups: structural models with conceptual correspondence to the AMM
  canonical form (mgcv, rstanarm, brms) and predictive agnostic models
  (ranger, xgboost). Both groups appear in the same script for
  convenience but should be interpreted under different lenses.
- The CRPS metric is computed via `scoringRules::crps_sample` when
  available, with a slower fallback implementation otherwise.
- Detailed methodological commentary is printed at the end of each
  benchmark section.

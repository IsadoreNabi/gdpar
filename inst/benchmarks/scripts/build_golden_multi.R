## build_golden_multi.R
##
## Generates the canonical multi-coord golden fits for Phase I.3 of
## Bloque 5.2: smoke_p2_auto, smoke_p2_cp, smoke_p2_ncp. Each golden
## is a single .rds saved to inst/benchmarks/golden/ following the
## four-layer schema produced by gdpar_snapshot_fit() and consumed by
## gdpar_golden_compare() (both defined in R/golden_compare.R):
##
##   structural : class signatures, slot shapes, column names.
##   discrete   : integer-valued sampler diagnostics (n_divergent,
##                treedepth_max, n_leapfrog_total).
##   continuous : posterior means + sd + Monte Carlo SE + ess_bulk
##                for theta_ref, a_coef, sigma_a, sigma_y, lp__.
##   sanity     : aggregate gates that must hold absolutely
##                (rhat_max, ess_bulk_min, ess_tail_min,
##                divergent_pct, ebfmi_min).
##
## Usage:
##   Rscript inst/benchmarks/scripts/build_golden_multi.R
##
## Wall time: about 5-7 minutes on the reference machine (one Stan
## compile multi shared across the three fits + three short fits
## with 1000/1000 iterations). The script delegates the snapshot
## extraction to gdpar_snapshot_fit() so the test suite and the
## generator stay in lockstep.

suppressPackageStartupMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())) # nolint
  library(devtools)
  load_all(".", quiet = TRUE)
})

build_smoke_p2_data <- function(seed = 42L, n = 300L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  beta_a_true <- list(c(2.0, -1.0), c(-2.0, 1.0))
  theta_ref_true <- c(0.5, -0.5)
  sigma_true <- 0.3
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1_c +
      beta_a_true[[1]][2] * x2_c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1_c +
      beta_a_true[[2]][2] * x2_c
  )
  y_mat <- eta_true + matrix(rnorm(n * 2L, 0, sigma_true), n, 2L)
  df <- data.frame(x1 = x1, x2 = x2)
  df$y <- y_mat
  df
}

build_golden <- function(scenario_id, parametrization, df, seed) {
  cat(sprintf("Fitting %s (parametrization=%s)...\n",
              scenario_id, parametrization))
  spec <- amm_spec(p = 2L,
                    dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  fit <- suppressMessages(suppressWarnings(
    gdpar(
      formula       = y ~ x1 + x2,
      family        = fam, amm = spec, data = df,
      chains        = 2L, iter_warmup = 1000L,
      iter_sampling = 1000L, refresh = 0L, seed = seed,
      parametrization = parametrization,
      skip_id_check = FALSE, verbose = FALSE
    )
  ))
  snapshot <- gdpar_snapshot_fit(fit)
  c(
    list(
      schema_version  = 1L,
      scenario_id     = scenario_id,
      parametrization = parametrization,
      gdpar_version   = as.character(utils::packageVersion("gdpar")),
      cmdstan_version = tryCatch(cmdstanr::cmdstan_version(),
                                  error = function(e) NA_character_),
      R_version       = R.version.string,
      seed            = seed,
      n               = nrow(df),
      p               = 2L,
      n_warmup        = 1000L,
      n_sampling      = 1000L,
      n_chains        = 2L,
      generated_at    = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    snapshot
  )
}

# --- main -------------------------------------------------------------

if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Package 'posterior' is required to build the goldens.")
}
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Package 'cmdstanr' is required to build the goldens.")
}

df <- build_smoke_p2_data(seed = 42L, n = 300L)
out_dir <- "inst/benchmarks/golden"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scenarios <- list(
  list(id = "smoke_p2_auto", parametrization = "auto"),
  list(id = "smoke_p2_cp",   parametrization = "cp"),
  list(id = "smoke_p2_ncp",  parametrization = "ncp")
)

for (s in scenarios) {
  g <- build_golden(s$id, s$parametrization, df, seed = 42L)
  out_path <- file.path(out_dir, paste0(s$id, ".rds"))
  saveRDS(g, out_path)
  cat(sprintf(
    "Wrote: %s  n_div=%d  rhat_max=%.4f  ess_bulk_min=%.0f  div_pct=%.4f\n",
    out_path, g$discrete$n_divergent, g$sanity$rhat_max,
    g$sanity$ess_bulk_min, g$sanity$divergent_pct
  ))
}

cat("Done.\n")

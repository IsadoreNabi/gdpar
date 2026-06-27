## Block RG, sub-phase RG.7 -- DRIVER that closes the four missing Tweedie cells
## of benchmark 9.2.O (N__N1__gdpar__{nef,mac,app,grl}) with the CERTIFIED-LIMIT
## Laplace ELPD (see inst/benchmarks/scripts/rg7_laplace_elpd.R for the verdict
## and the method). Each cell = coord-1 BETA (occurrence; fit normally by NUTS,
## the well-behaved coordinate) + coord-2 TWEEDIE (count; the near-deterministic
## canyon -> Laplace ELPD with a sensitivity check). Writes the cell RDS in the
## exact format the 9.2.O reconciler expects (obj$manifest_row), so a subsequent
## run of run_revalidation_9_2_O.R reconciles the 4 new rows into the manifest
## and the coverage gate reaches 80/80.
##
## Checkpointed/resumable per cell. Gated by GDPAR_RG7_LAPLACE=1.
##   GDPAR_RG7_LAPLACE_FORCE=1        recompute even if the cell RDS exists.
##   GDPAR_RG7_LAPLACE_SUBREGIONS     comma list of sub-regions (subset).
##   GDPAR_RG7_LAPLACE_S              Laplace draws (default 2000).
##   GDPAR_RG7_LAPLACE_SMOKE=1        cheap settings (small S, short refit).
##
## Usage (from the gdpar working directory; setsid + TMPDIR outside tmpfs):
##   GDPAR_RG7_LAPLACE=1 Rscript inst/benchmarks/run_rg7_laplace_elpd.R

if (Sys.getenv("GDPAR_RG7_LAPLACE") != "1") {
  message("GDPAR_RG7_LAPLACE is not 1; skipping. Set it to 1 to execute.")
  quit(save = "no", status = 0L)
}

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all(".", quiet = TRUE))
repo_root <- normalizePath(".")
bench_dir <- file.path(repo_root, "inst", "benchmarks")
scripts_dir <- file.path(bench_dir, "scripts")
source(file.path(bench_dir, "ebird_metrics_9_2_O.R"))
source(file.path(bench_dir, "scenarios_ebird_9_2_O.R"))
source(file.path(scripts_dir, "ebird_adapters_9_2_O.R"))
source(file.path(scripts_dir, "rg7_geometry_helpers.R"))
source(file.path(scripts_dir, "rg7_laplace_elpd.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a
results_dir <- file.path(bench_dir, "results")
cells_dir <- file.path(results_dir, "cells_9_2_O")
laplace_dir <- file.path(results_dir, "rg7_laplace")
log_path <- file.path(results_dir, "run_rg7_laplace_elpd.log")
dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(laplace_dir, showWarnings = FALSE, recursive = TRUE)

.elog <- function(...) {
  msg <- sprintf(...)
  cat(format(Sys.time(), "%H:%M:%S"), " ", msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n",
      sep = "", file = log_path, append = TRUE)
}

## --- manifest row builder (replicated EXACTLY from run_revalidation_9_2_O.R so
##     the reconciler appends an identical-schema row) ---------------------
.manifest_cols <- c(
  "cell_id", "phase", "track", "level", "method", "sub_region",
  "n_train", "n_test", "status", "wall_seconds",
  "y1_family", "y1_status", "y1_elpd_orig", "y1_elpd_loo", "y1_coverage90",
  "y1_pareto_k_max", "y1_nonfinite",
  "y2_family", "y2_status", "y2_elpd_orig", "y2_elpd_loo", "y2_coverage90",
  "y2_pareto_k_max", "y2_nonfinite",
  "rhat_max", "ess_bulk_min", "n_divergences", "error_class", "seed",
  "gdpar_version", "cmdstan_version", "timestamp")
.gdpar_version <- as.character(utils::packageVersion("gdpar"))
.cmdstan_version <- tryCatch(as.character(cmdstanr::cmdstan_version()),
                             error = function(e) NA_character_)
.manifest_row <- function(meta, cell) {
  c1 <- cell$coords[["1"]]; c2 <- cell$coords[["2"]]
  ec <- c(c1$error_class, c2$error_class); ec <- ec[!is.na(ec)]
  data.frame(
    cell_id = meta$cell_id, phase = "9.2.O", track = meta$track,
    level = meta$level, method = meta$method, sub_region = meta$sub_region,
    n_train = meta$n_train, n_test = meta$n_test, status = cell$status,
    wall_seconds = cell$wall_seconds,
    y1_family = c1$family_tag, y1_status = c1$status,
    y1_elpd_orig = c1$elpd %||% NA_real_, y1_elpd_loo = c1$elpd_loo %||% NA_real_,
    y1_coverage90 = c1$coverage90 %||% NA_real_,
    y1_pareto_k_max = c1$pareto_k_max %||% NA_real_,
    y1_nonfinite = c1$n_nonfinite %||% NA_integer_,
    y2_family = c2$family_tag, y2_status = c2$status,
    y2_elpd_orig = c2$elpd %||% NA_real_, y2_elpd_loo = c2$elpd_loo %||% NA_real_,
    y2_coverage90 = c2$coverage90 %||% NA_real_,
    y2_pareto_k_max = c2$pareto_k_max %||% NA_real_,
    y2_nonfinite = c2$n_nonfinite %||% NA_integer_,
    rhat_max = cell$diag$rhat_max %||% NA_real_,
    ess_bulk_min = cell$diag$ess_bulk_min %||% NA_real_,
    n_divergences = cell$diag$n_divergences %||% NA_real_,
    error_class = if (length(ec)) paste(ec, collapse = ";") else NA_character_,
    seed = meta$fit_seed, gdpar_version = .gdpar_version,
    cmdstan_version = .cmdstan_version,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE)
}

## --- settings -------------------------------------------------------
smoke <- Sys.getenv("GDPAR_RG7_LAPLACE_SMOKE") == "1"
force <- Sys.getenv("GDPAR_RG7_LAPLACE_FORCE") == "1"
S_draws <- {
  raw <- Sys.getenv("GDPAR_RG7_LAPLACE_S")
  if (nzchar(raw)) as.integer(raw) else if (smoke) 300L else 2000L
}
refit_iters <- if (smoke) {
  list(chains = 1L, iter_warmup = 150L, iter_sampling = 150L)
} else {
  list(chains = 2L, iter_warmup = 250L, iter_sampling = 300L)
}
refit_td <- if (smoke) 7L else 8L
subset_sr <- Sys.getenv("GDPAR_RG7_LAPLACE_SUBREGIONS")
sub_regions <- .ebird_9_2_O_sub_regions
if (nzchar(subset_sr)) {
  short2full <- setNames(names(.ebird_9_2_O_subregion_short),
                         unname(.ebird_9_2_O_subregion_short))
  req <- trimws(strsplit(subset_sr, ",")[[1L]])
  sub_regions <- unique(c(sub_regions[sub_regions %in% req],
                          short2full[req[req %in% names(short2full)]]))
}

.elog("=== RG7 LAPLACE BATCH START === S=%d smoke=%s subregions=%s",
      S_draws, smoke, paste(sub_regions, collapse = "/"))
ds <- ebird_9_2_O_data_cache(results_dir, scripts_dir, force = FALSE)

for (sr in sub_regions) {
  srs <- unname(.ebird_9_2_O_subregion_short[sr])
  cell_id <- sprintf("N__N1__gdpar__%s", srs)
  rds_path <- file.path(cells_dir, paste0(cell_id, ".rds"))
  if (!force && file.exists(rds_path)) {
    .elog("SKIP %s (exists)", cell_id); next
  }
  .elog("--- CELL %s START ---", cell_id)
  sp <- ebird_9_2_O_subregion_split(ds$full, sr)
  fit_seed <- ebird_9_2_O_fit_seed(sr, "N1", "gdpar")
  meta <- list(cell_id = cell_id, track = "N", level = "N1", method = "gdpar",
               sub_region = sr, fit_seed = fit_seed,
               n_train = nrow(sp$train), n_test = nrow(sp$test))

  ## cached mode (warm start) from the B9.32 geometry checkpoint, when valid.
  cached <- tryCatch(readRDS(file.path(results_dir, "rg7_geometry",
                                       sprintf("rg7_%s.rds", srs))),
                     error = function(e) NULL)
  cm <- tryCatch({
    v <- cached$mode$mode
    if (is.numeric(v) && length(v) == 14L) list(mode = as.numeric(v)) else NULL
  }, error = function(e) NULL)

  cell <- tryCatch({
    ## coord 1: beta occurrence (level Nb) -- fit normally by NUTS.
    .elog("%s coord1 beta (Nb) fit ...", cell_id)
    r1 <- .ebird_fit_gdpar_scalar(1L, "N", "Nb", sp$train, sp$test,
                                  fit_seed + 1L, "beta", "beta")
    c1 <- ebird_finalize_coord(r1, 1L, "N", sp$train, sp$test, fit_seed + 1L)

    ## coord 2: tweedie count (level Nt) -- Laplace ELPD.
    .elog("%s coord2 tweedie (Nt) Laplace ELPD ...", cell_id)
    lap <- rg7_laplace_tweedie_coord(sp$train, sp$test, fit_seed, cached = cm,
                                     S = S_draws, refit_iters = refit_iters,
                                     refit_treedepth = refit_td, verbose = TRUE)
    mp <- lap$mode_precision
    r2 <- list(status = "ok", family_tag = "tweedie",
               param_draws = lap$param_draws, loglik_native = NULL, loo = NULL,
               diag = list(rhat_max = NA_real_, ess_bulk_min = NA_real_,
                           n_divergences = NA_real_),
               wall_seconds = lap$wall_seconds, error_class = NA_character_)
    c2 <- ebird_finalize_coord(r2, 2L, "N", sp$train, sp$test, fit_seed + 2L)

    wall <- sum(c(c1$wall_seconds, c2$wall_seconds), na.rm = TRUE)
    ## cell-level diagnostics: the BETA coord is the genuine NUTS fit, so its
    ## rhat/ess/divergences populate the manifest (the tweedie coord is Laplace-
    ## certified, recorded honestly in the Laplace sidecar + the report).
    diag <- if (!is.null(c1$diag) && !all(is.na(unlist(c1$diag)))) c1$diag
            else .ebird_diag_null()
    cell <- .ebird_assemble_cell(c1, c2, wall, diag)
    attr(cell, "laplace") <- list(
      mode_grad = mp$grad_norm, mode_offset_sd = mp$mode_offset_sd,
      mode_converged = mp$converged, cond = mp$cond, eig_min = mp$eig_min,
      eig_max = mp$eig_max, all_pos = mp$all_pos, hess_method = mp$method,
      mode_source = mp$mode_source, verify_max_abs_diff = lap$verify_flatten$max_abs_diff,
      n_draw_ok = lap$n_draw_ok, n_draw_fail = lap$n_draw_fail,
      fail_rate = lap$fail_rate, S = S_draws,
      sensitivity = lap$sensitivity, fit_quality = lap$fit_quality,
      elpd_laplace = c2$elpd, elpd_plugin = lap$sensitivity$elpd_plugin,
      se_elpd = lap$sensitivity$se_elpd, coverage90 = c2$coverage90)
    cell
  }, error = function(e) {
    .elog("%s ERROR: %s", cell_id, conditionMessage(e))
    list(status = "fit_failed", wall_seconds = NA_real_,
         diag = .ebird_diag_null(),
         coords = list(
           `1` = list(status = "fit_failed", family_tag = "beta",
                      error_class = class(e)[1L]),
           `2` = list(status = "fit_failed", family_tag = "tweedie",
                      error_class = class(e)[1L])))
  })

  mrow <- .manifest_row(meta, cell)
  obj <- list(cell_meta = meta, manifest_row = mrow, cell = cell,
              status = cell$status, laplace_meta = attr(cell, "laplace"))
  tmp <- paste0(rds_path, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, rds_path)
  ## a stand-alone Laplace sidecar for the report (full certification record).
  saveRDS(obj$laplace_meta, file.path(laplace_dir, sprintf("laplace_%s.rds", srs)))
  .elog("=== CELL %s DONE === status=%s wall=%.0fs y1(beta)=%.2f[%s] y2(tweedie/laplace)=%.2f[%s]",
        cell_id, cell$status, cell$wall_seconds %||% NA_real_,
        cell$coords[["1"]]$elpd %||% NA_real_, cell$coords[["1"]]$status,
        cell$coords[["2"]]$elpd %||% NA_real_, cell$coords[["2"]]$status)
}
.elog("=== RG7 LAPLACE BATCH COMPLETE ===")

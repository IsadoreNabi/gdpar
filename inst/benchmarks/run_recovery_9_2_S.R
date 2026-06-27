## Block 9.2.S driver: external synthetic recovery benchmark vs the
## roster of 5 methods (gdpar + mgcv + brms + INLA + rstanarm), over the
## 64-scenario factorial x R repetitions. Checkpointed/resumable per
## cell (DESIGN_9_2_EXTERNAL section 3); the coverage gate (full grid
## with a recorded status per cell) governs closure, not any single run.
##
## Grid (Decision D, corrected B9.11 emergent D79): 64 scenarios x 5
## methods x R=10 = 3200 cells. The "4" in early drafts (2560) predated
## adding rstanarm to the roster; the count is roster-driven here.
##
## Env gates:
##   GDPAR_BENCH_9_2=1          required to run (else skips, exit 0).
##   GDPAR_BENCH_9_2_FORCE=1    recompute even if a cell RDS exists.
##   GDPAR_BENCH_9_2_SCENARIOS  comma list of scenario_ids (subset).
##   GDPAR_BENCH_9_2_METHODS    comma list of methods (subset).
##   GDPAR_BENCH_9_2_REPS       comma list of rep integers (subset).
##   GDPAR_BENCH_9_2_MAX_CELLS  stop after this many newly-computed cells
##                               (lote-por-sesion budget; default Inf).
##   GDPAR_BENCH_9_2_ITER       "chains,warmup,sampling" HMC budget
##                               (default 2,500,1000 = Block 7 setting).
##   GDPAR_BENCH_9_2_N          per-scenario n (default 2000).
##   GDPAR_BENCH_9_2_FAMILIES   comma list of gdpar_family (subset).
##   GDPAR_BENCH_9_2_RDS_ONLY=1 concurrency-safe mode (B9.12): write only
##                               the per-cell RDS to the canonical cells_dir
##                               and skip both the manifest reconcile and the
##                               manifest CSV append. Lets a second batch run
##                               concurrently with a first one on the SAME
##                               canonical cells_9_2_S/ (cell_ids never
##                               collide; RDS writes are atomic) without two
##                               processes appending to the same CSV. The
##                               manifest is regenerated from the RDS by the
##                               next non-RDS-only invocation's reconcile.
##
## Usage (from the gdpar working directory):
##   GDPAR_BENCH_9_2=1 Rscript inst/benchmarks/run_recovery_9_2_S.R

if (Sys.getenv("GDPAR_BENCH_9_2") != "1") {
  message("GDPAR_BENCH_9_2 is not 1; skipping the 9.2.S recovery run. ",
          "Set GDPAR_BENCH_9_2=1 to execute it.")
  quit(save = "no", status = 0L)
}

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to load gdpar via load_all().")
}
suppressMessages(devtools::load_all(".", quiet = TRUE))

repo_root <- normalizePath(".")
bench_dir <- file.path(repo_root, "inst", "benchmarks")
scripts_dir <- file.path(bench_dir, "scripts")
source(file.path(scripts_dir, "bench_multi_helpers.R"))
source(file.path(bench_dir, "recovery_metrics.R"))
source(file.path(bench_dir, "scenarios_recovery_factorial.R"))
source(file.path(scripts_dir, "recovery_adapters.R"))

results_dir <- file.path(bench_dir, "results")
## GDPAR_BENCH_9_2_TAG isolates a smoke/scratch run from the canonical
## artefacts (e.g. tag="smoke" -> cells_9_2_S_smoke/ + _smoke.csv).
.tag <- Sys.getenv("GDPAR_BENCH_9_2_TAG")
.suffix <- if (nzchar(.tag)) paste0("_", .tag) else ""
cells_dir   <- file.path(results_dir, paste0("cells_9_2_S", .suffix))
manifest_path <- file.path(results_dir, paste0("manifest_9_2_S", .suffix, ".csv"))
log_path <- file.path(results_dir, paste0("run_recovery_9_2_S", .suffix, ".log"))
dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)

.rec_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n",
      sep = "", file = log_path, append = TRUE)
}

## --- HMC budget / sizing -------------------------------------------
.parse_iter <- function() {
  raw <- Sys.getenv("GDPAR_BENCH_9_2_ITER")
  if (!nzchar(raw)) return(list(chains = 2L, iter_warmup = 500L,
                                iter_sampling = 1000L))
  v <- as.integer(strsplit(raw, ",", fixed = TRUE)[[1L]])
  list(chains = v[1L], iter_warmup = v[2L], iter_sampling = v[3L])
}
iter <- .parse_iter()
n_total <- {
  raw <- Sys.getenv("GDPAR_BENCH_9_2_N")
  if (nzchar(raw)) as.integer(raw) else 2000L
}
R <- 10L

## --- manifest schema ------------------------------------------------
.manifest_cols <- c(
  "cell_id", "phase", "regime", "bitmask", "pathologies_active",
  "method", "rep", "seed", "gdpar_family", "K", "n", "status",
  "wall_seconds", "rmse_mean", "bias_mean", "rmse_sd", "tail_metric",
  "coverage90", "elpd_loo", "pareto_k_max", "rhat_max", "ess_bulk_min",
  "n_divergences", "error_class", "gdpar_version", "cmdstan_version",
  "timestamp")

.gdpar_version <- as.character(utils::packageVersion("gdpar"))
.cmdstan_version <- tryCatch(
  as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_)

.manifest_row <- function(cell_meta, mrow) {
  data.frame(
    cell_id = cell_meta$cell_id, phase = "9.2.S",
    regime = cell_meta$regime, bitmask = cell_meta$bitmask,
    pathologies_active = cell_meta$pathologies_active,
    method = cell_meta$method, rep = cell_meta$rep,
    seed = cell_meta$fit_seed, gdpar_family = cell_meta$gdpar_family,
    K = cell_meta$K, n = cell_meta$n, status = mrow$status,
    wall_seconds = mrow$wall_seconds, rmse_mean = mrow$rmse_mean,
    bias_mean = mrow$bias_mean, rmse_sd = mrow$rmse_sd,
    tail_metric = mrow$tail_metric, coverage90 = mrow$coverage90,
    elpd_loo = mrow$elpd_loo, pareto_k_max = mrow$pareto_k_max,
    rhat_max = mrow$rhat_max, ess_bulk_min = mrow$ess_bulk_min,
    n_divergences = mrow$n_divergences, error_class = mrow$error_class,
    gdpar_version = .gdpar_version, cmdstan_version = .cmdstan_version,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE)
}

.append_manifest <- function(row) {
  write_header <- !file.exists(manifest_path)
  utils::write.table(row[, .manifest_cols], file = manifest_path,
                     sep = ",", row.names = FALSE,
                     col.names = write_header, append = !write_header,
                     qmethod = "double")
}

## Reconcile manifest with RDS on disk (defend against a crash between
## the atomic rename and the manifest append).
.reconcile_manifest <- function() {
  rds_files <- list.files(cells_dir, pattern = "\\.rds$", full.names = TRUE)
  done_ids <- sub("\\.rds$", "", basename(rds_files))
  recorded <- character(0)
  if (file.exists(manifest_path)) {
    mm <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                          colClasses = "character")
    recorded <- mm$cell_id
  }
  missing <- setdiff(done_ids, recorded)
  for (id in missing) {
    obj <- tryCatch(readRDS(file.path(cells_dir, paste0(id, ".rds"))),
                    error = function(e) NULL)
    if (!is.null(obj) && !is.null(obj$manifest_row)) {
      .append_manifest(obj$manifest_row)
    }
  }
  if (length(missing) > 0L) {
    .rec_log("Reconciled %d RDS cells missing from manifest.",
             length(missing))
  }
  invisible(NULL)
}

## --- cell grid ------------------------------------------------------
tab <- recovery_scenario_table()
methods <- names(.rec_method_idx)
sub_sc  <- Sys.getenv("GDPAR_BENCH_9_2_SCENARIOS")
sub_me  <- Sys.getenv("GDPAR_BENCH_9_2_METHODS")
sub_rep <- Sys.getenv("GDPAR_BENCH_9_2_REPS")
sub_reg <- Sys.getenv("GDPAR_BENCH_9_2_REGIME")   ## "C" or "K"
sub_fam <- Sys.getenv("GDPAR_BENCH_9_2_FAMILIES") ## comma list of gdpar_family
if (nzchar(sub_reg)) {
  tab <- tab[tab$regime == sub_reg, , drop = FALSE]
}
if (nzchar(sub_fam)) {
  fams <- trimws(strsplit(sub_fam, ",", fixed = TRUE)[[1L]])
  tab <- tab[tab$gdpar_family %in% fams, , drop = FALSE]
}
if (nzchar(sub_sc)) {
  want <- trimws(strsplit(sub_sc, ",", fixed = TRUE)[[1L]])
  tab <- tab[tab$scenario_id %in% want, , drop = FALSE]
}
if (nzchar(sub_me)) {
  methods <- intersect(methods,
                       trimws(strsplit(sub_me, ",", fixed = TRUE)[[1L]]))
}
reps <- if (nzchar(sub_rep)) {
  as.integer(trimws(strsplit(sub_rep, ",", fixed = TRUE)[[1L]]))
} else seq_len(R)

max_cells <- {
  raw <- Sys.getenv("GDPAR_BENCH_9_2_MAX_CELLS")
  if (nzchar(raw)) as.integer(raw) else .Machine$integer.max
}
force <- Sys.getenv("GDPAR_BENCH_9_2_FORCE") == "1"
rds_only <- Sys.getenv("GDPAR_BENCH_9_2_RDS_ONLY") == "1"

method_fns <- recovery_method_table()

if (!rds_only) .reconcile_manifest()

.rec_log("=== 9_2_S BATCH START === scenarios=%d methods=%s reps=%s n=%d iter=%s rds_only=%s",
         nrow(tab), paste(methods, collapse = "/"),
         paste(reps, collapse = "/"), n_total,
         paste(unlist(iter), collapse = ","), rds_only)

computed <- 0L; skipped <- 0L
for (si in seq_len(nrow(tab))) {
  srow <- tab[si, , drop = FALSE]
  for (rep in reps) {
    realized <- NULL                 ## lazily realize per (scenario, rep)
    for (m in methods) {
      cell_id <- sprintf("%s__%s__rep%02d", srow$scenario_id, m, rep)
      rds_path <- file.path(cells_dir, paste0(cell_id, ".rds"))
      if (!force && file.exists(rds_path)) { skipped <- skipped + 1L; next }
      if (computed >= max_cells) {
        .rec_log("=== 9_2_S BATCH BUDGET REACHED === computed=%d", computed)
        .rec_log("=== 9_2_S BATCH COMPLETE === computed=%d skipped=%d",
                 computed, skipped)
        quit(save = "no", status = 0L)
      }
      if (is.null(realized)) {
        realized <- recovery_realize_cell(srow, rep, n = n_total)
      }
      fit_seed <- rec_fit_seed(srow$truth_seed, rep, .rec_method_idx[[m]])
      cell_meta <- list(
        cell_id = cell_id, scenario_id = srow$scenario_id,
        regime = srow$regime, bitmask = srow$bitmask,
        pathologies_active = srow$pathologies_active, method = m, rep = rep,
        fit_seed = fit_seed, gdpar_family = srow$gdpar_family,
        K = srow$K, n = n_total)
      pred <- tryCatch(
        method_fns[[m]](realized$model_plan, realized$train, realized$test,
                        realized$grid, fit_seed, iter),
        error = function(e) rec_pred_status(
          "fit_failed", error_class = class(e)[1L],
          gdpar_family = srow$gdpar_family, K = srow$K))
      mrow <- recovery_compute_metrics(pred, realized$truth_grid)
      manifest_row <- .manifest_row(cell_meta, mrow)
      obj <- list(cell_meta = cell_meta, metrics = mrow,
                  manifest_row = manifest_row,
                  diag = pred$diag, status = pred$status,
                  truth_summary = list(
                    mu_mean = mean(realized$truth_grid$mu),
                    sd_mean = mean(realized$truth_grid$sd)))
      tmp <- paste0(rds_path, ".tmp")
      saveRDS(obj, tmp); file.rename(tmp, rds_path)
      if (!rds_only) .append_manifest(manifest_row)
      computed <- computed + 1L
      .rec_log("=== 9_2_S CELL %s DONE === status=%s wall=%.1fs",
               cell_id, mrow$status, mrow$wall_seconds %||% NA_real_)
    }
    rm(realized); gc(verbose = FALSE)
  }
}

.rec_log("=== 9_2_S BATCH COMPLETE === computed=%d skipped=%d",
         computed, skipped)

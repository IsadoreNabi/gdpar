## Block RG, RG.1.c driver: full calibration of the posterior-geometry
## diagnostic classifier against gdpar_geometry_suite. Replaces the "8/8
## diagonal of one replica" smoke (B9.22) with measured error rates over a
## (target x difficulty x pilot-budget x replica) grid. Checkpointed/resumable
## per cell (the 9.2.S pattern); coverage is governed by a recorded status per
## cell, not by any single run.
##
## The expensive part (collecting signals via cheap NUTS pilots) is done here;
## the threshold sweep is a cheap offline pass over the stored signals
## (rg_calibration_analysis.R), so this driver only collects and never tunes.
##
## Env gates:
##   GDPAR_RG=1              required to run (else skips, exit 0).
##   GDPAR_RG_FORCE=1        recompute even if a cell RDS exists.
##   GDPAR_RG_TAG=...        isolate a smoke/scratch run (tag="smoke" ->
##                            cells_rg_smoke/ + manifest_rg_smoke.csv).
##   GDPAR_RG_REPS=...       comma list or "a:b" range of replica indices
##                            (default 1:5 = Stage 1). Stage 2 passes the
##                            top-up reps for the boundary cells.
##   GDPAR_RG_CELLS=...      comma list of cell_base values (subset; default
##                            all 48). Stage 2 passes the boundary cells.
##   GDPAR_RG_MAX_CELLS=n    stop after n newly-computed cells (lote budget).
##
## Usage (from the gdpar working directory):
##   GDPAR_RG=1 Rscript inst/benchmarks/run_rg_calibration.R

if (Sys.getenv("GDPAR_RG") != "1") {
  message("GDPAR_RG is not 1; skipping the RG.1.c calibration run. ",
          "Set GDPAR_RG=1 to execute it.")
  quit(save = "no", status = 0L)
}

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to load gdpar via load_all().")
}
suppressMessages(devtools::load_all(".", quiet = TRUE))

repo_root   <- normalizePath(".")
bench_dir   <- file.path(repo_root, "inst", "benchmarks")
source(file.path(bench_dir, "scenarios_geometry_calibration.R"))

results_dir <- file.path(bench_dir, "results")
.tag    <- Sys.getenv("GDPAR_RG_TAG")
.suffix <- if (nzchar(.tag)) paste0("_", .tag) else ""
cells_dir     <- file.path(results_dir, paste0("cells_rg", .suffix))
manifest_path <- file.path(results_dir, paste0("manifest_rg", .suffix, ".csv"))
log_path      <- file.path(results_dir, paste0("run_rg_calibration", .suffix, ".log"))
dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)

.rg_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n",
      sep = "", file = log_path, append = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.gdpar_version  <- as.character(utils::packageVersion("gdpar"))
.cmdstan_version <- tryCatch(
  as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_)

## --- compact manifest (full signals live in the per-cell RDS) -------
.manifest_cols <- c(
  "cell_id", "cell_base", "target_id", "diff_level", "diff_value",
  "pilot_level", "warmup", "sampling", "rep", "seed", "ground_truth",
  "default_pathology", "default_correct", "slope", "status", "elapsed",
  "gdpar_version", "cmdstan_version", "timestamp")

.manifest_row <- function(obj) {
  data.frame(
    cell_id = obj$cell_id, cell_base = obj$cell_base,
    target_id = obj$target_id, diff_level = obj$diff_level,
    diff_value = obj$diff_value, pilot_level = obj$pilot_level,
    warmup = obj$warmup, sampling = obj$sampling, rep = obj$rep,
    seed = obj$seed, ground_truth = obj$ground_truth,
    default_pathology = obj$default_pathology,
    default_correct = obj$default_correct,
    slope = obj$slope %||% NA_real_, status = obj$status,
    elapsed = obj$elapsed,
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

## Reconcile manifest with RDS on disk (defend against a crash between the
## atomic rename and the manifest append).
.reconcile_manifest <- function() {
  rds_files <- list.files(cells_dir, pattern = "\\.rds$", full.names = TRUE)
  done_ids  <- sub("\\.rds$", "", basename(rds_files))
  recorded  <- character(0)
  if (file.exists(manifest_path)) {
    mm <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                          colClasses = "character")
    recorded <- mm$cell_id
  }
  missing <- setdiff(done_ids, recorded)
  for (id in missing) {
    obj <- tryCatch(readRDS(file.path(cells_dir, paste0(id, ".rds"))),
                    error = function(e) NULL)
    if (!is.null(obj)) .append_manifest(.manifest_row(obj))
  }
  if (length(missing) > 0L) {
    .rg_log("Reconciled %d RDS cells missing from manifest.", length(missing))
  }
  invisible(NULL)
}

## --- rep range / cell subset / budget -------------------------------
.parse_reps <- function() {
  raw <- Sys.getenv("GDPAR_RG_REPS")
  if (!nzchar(raw)) return(1:5)
  if (grepl(":", raw, fixed = TRUE)) {
    ab <- as.integer(strsplit(raw, ":", fixed = TRUE)[[1L]])
    return(seq.int(ab[1L], ab[2L]))
  }
  as.integer(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]]))
}
reps <- .parse_reps()

grid <- rg_cell_grid()
sub_cells <- Sys.getenv("GDPAR_RG_CELLS")
if (nzchar(sub_cells)) {
  want <- trimws(strsplit(sub_cells, ",", fixed = TRUE)[[1L]])
  grid <- grid[grid$cell_base %in% want, , drop = FALSE]
}

## Adaptive Stage-2 selection (GDPAR_RG_ADAPTIVE=1): restrict the grid to the
## cell_bases whose Wilson CI for the default-threshold correct rate is still
## UNRESOLVED -- it overlaps 0.5 (the decision boundary), or the rate is 0 (an
## honest-failure case whose precise rate matters for the report). This is the
## minimax-robust top-up: replicas are spent only where the verdict is
## statistically undecided; the already-settled cells (5/5 correct, or a clear
## miss bounded away from 0.5) keep their Stage-1 budget. Reproducible: computed
## from the existing RDS, not by hand.
if (Sys.getenv("GDPAR_RG_ADAPTIVE") == "1") {
  source(file.path(bench_dir, "rg_calibration_analysis.R"))
  tab_so_far <- rg_collect(cells_dir)
  th0 <- gdpar_geometry_thresholds()
  pred0 <- rg_predict(tab_so_far, th0)
  ok0 <- pred0 == tab_so_far$ground_truth
  by_base <- tapply(seq_along(ok0), tab_so_far$cell_base, function(idx) {
    x <- sum(ok0[idx]); n <- length(idx)
    w <- rg_wilson(x, n)
    rate <- x / n
    list(rate = rate, lo = w["lo"], hi = w["hi"])
  })
  keep <- names(by_base)[vapply(by_base, function(z) {
    overlaps_half <- isTRUE(z$lo <= 0.5 && z$hi >= 0.5)
    overlaps_half || isTRUE(z$rate == 0)
  }, logical(1))]
  grid <- grid[grid$cell_base %in% keep, , drop = FALSE]
  cat(sprintf("Adaptive Stage 2: %d unresolved cell_bases selected.\n",
              length(keep)))
}

max_cells <- {
  raw <- Sys.getenv("GDPAR_RG_MAX_CELLS")
  if (nzchar(raw)) as.integer(raw) else .Machine$integer.max
}
force <- Sys.getenv("GDPAR_RG_FORCE") == "1"

suite <- gdpar_geometry_suite()

.reconcile_manifest()
.rg_log("=== RG.1.c BATCH START === cell_bases=%d reps=%s (%d cells max=%s)",
        nrow(grid), paste(range(reps), collapse = ":"),
        nrow(grid) * length(reps),
        if (is.finite(max_cells)) max_cells else "Inf")

computed <- 0L; skipped <- 0L
for (ci in seq_len(nrow(grid))) {
  crow <- grid[ci, , drop = FALSE]
  for (rep in reps) {
    cell_id  <- sprintf("%s__rep%02d", crow$cell_base, rep)
    rds_path <- file.path(cells_dir, paste0(cell_id, ".rds"))
    if (!force && file.exists(rds_path)) { skipped <- skipped + 1L; next }
    if (computed >= max_cells) {
      .rg_log("=== RG.1.c BUDGET REACHED === computed=%d", computed)
      break
    }
    seed <- rg_cell_seed(crow$cell_base, rep)
    obj <- tryCatch(
      rg_run_cell(suite, crow, rep, seed),
      error = function(e) {
        .rg_log("CELL %s ERROR: %s", cell_id, conditionMessage(e))
        list(cell_id = cell_id, cell_base = crow$cell_base,
             target_id = crow$target_id, diff_level = crow$diff_level,
             diff_value = crow$diff_value, pilot_level = crow$pilot_level,
             warmup = crow$warmup, sampling = crow$sampling, rep = rep,
             seed = seed, signals = NULL, slope = NA_real_, culprit = NULL,
             ground_truth = suite[[crow$target_id]]$pathology,
             default_pathology = NA_character_, default_correct = NA,
             elapsed = NA_real_, status = "cell_error")
      })
    tmp <- paste0(rds_path, ".tmp")
    saveRDS(obj, tmp); file.rename(tmp, rds_path)
    .append_manifest(.manifest_row(obj))
    computed <- computed + 1L
    .rg_log("=== RG.1.c CELL %s DONE === gt=%s pred=%s correct=%s wall=%.1fs",
            cell_id, obj$ground_truth, obj$default_pathology,
            obj$default_correct, obj$elapsed %||% NA_real_)
  }
  if (computed >= max_cells) break
}

.rg_log("=== RG.1.c BATCH COMPLETE === computed=%d skipped=%d", computed, skipped)

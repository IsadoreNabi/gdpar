## Block 9.2.O driver: eBird external re-validation vs the roster of 5
## methods (gdpar + mgcv + brms + INLA + rstanarm), over the two encoding
## tracks (G gaussian-faithful K=1,2,3; N natural-support beta/tweedie) x
## 4 NE-USA sub-regions (DESIGN_9_2_EXTERNAL section 4, refined B9.19).
## Checkpointed/resumable per cell (same mechanism as 9.2.S); the coverage
## gate (every cell carrying a recorded status) governs closure.
##
## Grid: 80 cells = (G: 3 levels + N: 1 level) x 5 methods x 4 sub-regions.
## Real data -> ONE deterministic fit per cell, no reps.
##
## Env gates:
##   GDPAR_BENCH_9_2_O=1           required to run (else skips, exit 0).
##   GDPAR_BENCH_9_2_O_FORCE=1     recompute even if a cell RDS exists.
##   GDPAR_BENCH_9_2_O_TRACKS      comma list in {G,N} (subset).
##   GDPAR_BENCH_9_2_O_LEVELS      comma list in {G1,G2,G3,N1} (subset).
##   GDPAR_BENCH_9_2_O_METHODS     comma list of methods (subset).
##   GDPAR_BENCH_9_2_O_SUBREGIONS  comma list of sub-regions (subset).
##   GDPAR_BENCH_9_2_O_MAX_CELLS   stop after this many new cells (budget).
##   GDPAR_BENCH_9_2_O_RDS_ONLY=1  concurrency-safe (write only RDS, skip
##                                  manifest reconcile + append).
##   GDPAR_BENCH_9_2_O_TAG         isolate a smoke/scratch run.
##   GDPAR_BENCH_9_2_O_DATA_FORCE=1 rebuild the cached eBird dataset.
##
## Usage (from the gdpar working directory):
##   GDPAR_BENCH_9_2_O=1 Rscript inst/benchmarks/run_revalidation_9_2_O.R

if (Sys.getenv("GDPAR_BENCH_9_2_O") != "1") {
  message("GDPAR_BENCH_9_2_O is not 1; skipping the 9.2.O eBird run. ",
          "Set GDPAR_BENCH_9_2_O=1 to execute it.")
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
source(file.path(bench_dir, "ebird_metrics_9_2_O.R"))
source(file.path(bench_dir, "scenarios_ebird_9_2_O.R"))
source(file.path(scripts_dir, "ebird_adapters_9_2_O.R"))

results_dir <- file.path(bench_dir, "results")
.tag <- Sys.getenv("GDPAR_BENCH_9_2_O_TAG")
.suffix <- if (nzchar(.tag)) paste0("_", .tag) else ""
cells_dir <- file.path(results_dir, paste0("cells_9_2_O", .suffix))
manifest_path <- file.path(results_dir, paste0("manifest_9_2_O", .suffix, ".csv"))
log_path <- file.path(results_dir, paste0("run_revalidation_9_2_O", .suffix, ".log"))
dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)

.elog <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n",
      sep = "", file = log_path, append = TRUE)
}

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
  ec <- c(c1$error_class, c2$error_class)
  ec <- ec[!is.na(ec)]
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

.append_manifest <- function(row) {
  write_header <- !file.exists(manifest_path)
  utils::write.table(row[, .manifest_cols], file = manifest_path, sep = ",",
                     row.names = FALSE, col.names = write_header,
                     append = !write_header, qmethod = "double")
}

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
    if (!is.null(obj) && !is.null(obj$manifest_row)) .append_manifest(obj$manifest_row)
  }
  if (length(missing) > 0L) .elog("Reconciled %d RDS cells into manifest.",
                                  length(missing))
  invisible(NULL)
}

## --- data + grid ----------------------------------------------------
data_force <- Sys.getenv("GDPAR_BENCH_9_2_O_DATA_FORCE") == "1"
.elog("=== 9_2_O loading eBird data cache (force=%s) ===", data_force)
ds <- ebird_9_2_O_data_cache(results_dir, scripts_dir, force = data_force)
.elog("eBird data: %d rows, zero-share count = %.3f, sub-regions = %s",
      nrow(ds$full), ds$zero_share_count,
      paste(sort(unique(ds$full$sub_region)), collapse = "/"))

tab <- ebird_9_2_O_cell_table()
sub_tr <- Sys.getenv("GDPAR_BENCH_9_2_O_TRACKS")
sub_lv <- Sys.getenv("GDPAR_BENCH_9_2_O_LEVELS")
sub_me <- Sys.getenv("GDPAR_BENCH_9_2_O_METHODS")
sub_sr <- Sys.getenv("GDPAR_BENCH_9_2_O_SUBREGIONS")
if (nzchar(sub_tr)) tab <- tab[tab$track %in% trimws(strsplit(sub_tr, ",")[[1L]]), ]
if (nzchar(sub_lv)) tab <- tab[tab$level %in% trimws(strsplit(sub_lv, ",")[[1L]]), ]
if (nzchar(sub_me)) tab <- tab[tab$method %in% trimws(strsplit(sub_me, ",")[[1L]]), ]
if (nzchar(sub_sr)) tab <- tab[tab$sub_region %in% trimws(strsplit(sub_sr, ",")[[1L]]), ]

max_cells <- {
  raw <- Sys.getenv("GDPAR_BENCH_9_2_O_MAX_CELLS")
  if (nzchar(raw)) as.integer(raw) else .Machine$integer.max
}
force <- Sys.getenv("GDPAR_BENCH_9_2_O_FORCE") == "1"
rds_only <- Sys.getenv("GDPAR_BENCH_9_2_O_RDS_ONLY") == "1"

if (!rds_only) .reconcile_manifest()

.elog("=== 9_2_O BATCH START === cells=%d methods=%s rds_only=%s",
      nrow(tab), paste(sort(unique(tab$method)), collapse = "/"), rds_only)

## memoize per-sub-region splits (avoid re-splitting per cell)
.splits <- new.env(parent = emptyenv())
.get_split <- function(sr) {
  key <- sr
  if (is.null(.splits[[key]])) {
    .splits[[key]] <- ebird_9_2_O_subregion_split(ds$full, sr)
  }
  .splits[[key]]
}

computed <- 0L; skipped <- 0L
for (i in seq_len(nrow(tab))) {
  row <- tab[i, , drop = FALSE]
  cell_id <- row$cell_id
  rds_path <- file.path(cells_dir, paste0(cell_id, ".rds"))
  if (!force && file.exists(rds_path)) { skipped <- skipped + 1L; next }
  if (computed >= max_cells) {
    .elog("=== 9_2_O BATCH BUDGET REACHED === computed=%d", computed)
    break
  }
  sp <- .get_split(row$sub_region)
  meta <- list(cell_id = cell_id, track = row$track, level = row$level,
               method = row$method, sub_region = row$sub_region,
               fit_seed = row$fit_seed, n_train = nrow(sp$train),
               n_test = nrow(sp$test))
  cell <- tryCatch(
    ebird_run_cell_9_2_O(row$track, row$level, row$method, sp$train, sp$test,
                         row$fit_seed),
    error = function(e) list(
      status = "fit_failed", wall_seconds = NA_real_, diag = list(),
      coords = list(`1` = list(status = "fit_failed", family_tag = NA_character_,
                               error_class = class(e)[1L]),
                    `2` = list(status = "fit_failed", family_tag = NA_character_,
                               error_class = class(e)[1L]))))
  mrow <- .manifest_row(meta, cell)
  obj <- list(cell_meta = meta, manifest_row = mrow, cell = cell,
              status = cell$status)
  tmp <- paste0(rds_path, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, rds_path)
  if (!rds_only) .append_manifest(mrow)
  computed <- computed + 1L
  .elog("=== 9_2_O CELL %s DONE === status=%s wall=%.1fs y1=%s y2=%s",
        cell_id, cell$status, cell$wall_seconds %||% NA_real_,
        cell$coords[["1"]]$status, cell$coords[["2"]]$status)
}

.elog("=== 9_2_O BATCH COMPLETE === computed=%d skipped=%d", computed, skipped)

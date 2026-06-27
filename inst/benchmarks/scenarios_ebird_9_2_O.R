## Block 9.2.O: scenario/grid generator for the eBird external
## re-validation (DESIGN_9_2_EXTERNAL section 4, refined B9.19).
##
## Inherits the Block 7 eBird Status & Trends p = 2 dataset (5 NE-USA
## forest passerines x 4 sub-regions; logit-occurrence + log1p-count,
## three WorldClim predictors, state-level grouping) via
## build_ebird_path_a_multi_dataset() in scripts/ebird_loader.R. The full
## sampled table is built ONCE and cached to RDS (the rasters are ~2 GB);
## every cell reads the cache, so the run is resumable without reloading.
##
## Two encoding tracks (B9.19, user decision "implement both + max-
## robustness countermeasures"):
##   track G (gaussian-faithful): model z1 = logit(occ), z2 = log1p(cnt)
##     with the K = 1,2,3 ladder gaussian / gaussian-heteroscedastic /
##     student_t. Full roster; ELPD comparable on the original measure via
##     the Jacobian (ebird_metrics_9_2_O.R).
##   track N (natural support): model occurrence with beta (mu ~ x, phi
##     population) and count with TWEEDIE (mu ~ x; zero-robust compound
##     Poisson-gamma -- the count is continuous abundance with up to ~70%
##     zeros, so gamma is invalid; tweedie is the correct support and is
##     gdpar's own family, Block 8.3.5b). beta: full roster; tweedie:
##     gdpar + mgcv only (the only roster members with a native tweedie),
##     the rest recorded unsupported (coverage-gate pass).
##
## Structure (B9.19 Q2 = "nativa-mixta"): gdpar/brms fit p = 2 jointly,
## mgcv/INLA/rstanarm per coordinate -- exactly as Block 7. Fits are PER
## SUB-REGION (B9.19 Q3) so the per-sub-region ELPD needed for criterion
## 2.O comes directly from the cell, no decomposition.
##
## Cell = (track, level, method, sub_region). Grid:
##   G: 3 levels x 5 methods x 4 sub-regions = 60
##   N: 1 level  x 5 methods x 4 sub-regions = 20
##   total 80 cells. Real data -> ONE deterministic fit per cell, no reps.

.ebird_9_2_O_sub_regions <- c("NE_FOREST", "MID_ATL_COAST",
                              "APPALACHIAN", "GREAT_LAKES")
.ebird_9_2_O_subregion_short <- c(NE_FOREST = "nef", MID_ATL_COAST = "mac",
                                  APPALACHIAN = "app", GREAT_LAKES = "grl")
.ebird_9_2_O_methods <- c("gdpar", "mgcv", "brms", "inla", "rstanarm")
.ebird_9_2_O_covariates <- c("lon", "lat", "bio1", "bio12", "bio15")
.ebird_9_2_O_seed_base <- 91069L          ## reserved 91069-91099 (B9.19)

## level table: track, level id, per-coordinate family tag + K + a short
## human label. For track N the two coordinates use DIFFERENT families
## (beta on y1, tweedie on y2), so family is recorded per coordinate.
ebird_9_2_O_level_meta <- function(track, level) {
  if (track == "G") {
    switch(level,
      G1 = list(track = "G", level = "G1", K = 1L,
                fam1 = "gaussian", fam2 = "gaussian",
                tag1 = "gaussian_z", tag2 = "gaussian_z",
                label = "gaussian K1"),
      G2 = list(track = "G", level = "G2", K = 2L,
                fam1 = "gaussian", fam2 = "gaussian",
                tag1 = "gaussian_z", tag2 = "gaussian_z",
                label = "gaussian-hetero K2"),
      G3 = list(track = "G", level = "G3", K = 3L,
                fam1 = "student_t", fam2 = "student_t",
                tag1 = "student_t_z", tag2 = "student_t_z",
                label = "student_t K3"),
      stop("bad G level"))
  } else {
    list(track = "N", level = "N1", K = 2L,
         fam1 = "beta", fam2 = "tweedie",
         tag1 = "beta", tag2 = "tweedie",
         label = "beta/tweedie natural")
  }
}

## Build (or load) the cached full sampled eBird dataset. Combines the
## Block 7 loader's train+test into one table and adds the modelling
## columns used by both tracks: z1/z2 (track G) are the loader's y1/y2;
## occ/cnt (track N) are the raw original-scale outcomes.
ebird_9_2_O_data_cache <- function(results_dir, scripts_dir, force = FALSE) {
  cache_path <- file.path(results_dir, "ebird_9_2_O_data.rds")
  if (!force && file.exists(cache_path)) return(readRDS(cache_path))
  ## build_ebird_path_a_multi_dataset() (ebird_loader.R) depends on
  ## .family_spec_gaussian_homog() from synthetic_hard_multi.R -- the same
  ## dependency run_ebird_a_multi.R sources (Block 7). Source both here.
  bench_dir <- dirname(scripts_dir)
  source(file.path(bench_dir, "synthetic_hard_multi.R"), local = TRUE)
  source(file.path(scripts_dir, "ebird_loader.R"), local = TRUE)
  ds <- build_ebird_path_a_multi_dataset()
  full <- rbind(ds$train, ds$test)
  rownames(full) <- NULL
  full$z1  <- full$y1                      ## logit-occurrence (track G)
  full$z2  <- full$y2                      ## log1p-count      (track G)
  full$occ <- full$occurrence_mean         ## original [0,1]   (track N y1)
  full$cnt <- full$count_mean              ## original [0,Inf) (track N y2)
  full$group <- droplevels(factor(full$group))
  obj <- list(full = full, family_spec = ds$family_spec,
              strata_summary = ds$strata_summary,
              zero_share_count = mean(full$cnt == 0))
  tmp <- paste0(cache_path, ".tmp"); saveRDS(obj, tmp)
  file.rename(tmp, cache_path)
  obj
}

## Deterministic 80/20 split of one sub-region's rows.
ebird_9_2_O_subregion_split <- function(full, sub_region, train_frac = 0.8) {
  d <- full[full$sub_region == sub_region, , drop = FALSE]
  d$group <- droplevels(factor(d$group))
  sr_idx <- match(sub_region, .ebird_9_2_O_sub_regions)
  set.seed(.ebird_9_2_O_seed_base + sr_idx)
  n <- nrow(d); sh <- sample.int(n)
  d <- d[sh, , drop = FALSE]; rownames(d) <- NULL
  n_tr <- floor(train_frac * n)
  list(train = d[seq_len(n_tr), , drop = FALSE],
       test  = d[(n_tr + 1L):n, , drop = FALSE],
       sub_region = sub_region)
}

## Deterministic per-cell fit seed.
ebird_9_2_O_fit_seed <- function(sub_region, level, method) {
  sr  <- match(sub_region, .ebird_9_2_O_sub_regions)
  lv  <- match(level, c("G1", "G2", "G3", "N1"))
  me  <- match(method, .ebird_9_2_O_methods)
  .ebird_9_2_O_seed_base + 1000L * sr + 10L * lv + me
}

## The 80-cell grid.
ebird_9_2_O_cell_table <- function() {
  rows <- list()
  levels_by_track <- list(G = c("G1", "G2", "G3"), N = c("N1"))
  for (track in c("G", "N")) {
    for (level in levels_by_track[[track]]) {
      for (method in .ebird_9_2_O_methods) {
        for (sr in .ebird_9_2_O_sub_regions) {
          srs <- unname(.ebird_9_2_O_subregion_short[sr])
          cell_id <- sprintf("%s__%s__%s__%s", track, level, method, srs)
          rows[[length(rows) + 1L]] <- data.frame(
            cell_id = cell_id, track = track, level = level,
            method = method, sub_region = sr,
            fit_seed = ebird_9_2_O_fit_seed(sr, level, method),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  do.call(rbind, rows)
}

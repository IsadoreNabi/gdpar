## eBird Status & Trends loader for the Block 7 external-validation
## bench. Builds a p = 2 multivariate dataset from cell-wise seasonal
## (breeding) outcomes of five NE-USA forest passerines, joined with
## three WorldClim 2.1 bioclimatic predictors and a state-level
## grouping factor used both as gdpar grouping and as a random-effect
## term in the four competitor wrappers.
##
## Exposes:
##   .NE_USA_SUBREGION_MAP            -> character vector state -> region
##   .assign_state_to_subregion(s)    -> region label or NA
##   .load_ebird_species_seasonal_p2(sp, states_sf, worldclim_rast,
##                                    resolution = "9km")
##                                    -> data.frame
##   .grid_sample_stratified(df, strata_col, n_max)
##                                    -> data.frame
##   build_ebird_path_a_multi_dataset(species = ...,
##                                    resolution = "9km",
##                                    n_max_per_strata = 750L,
##                                    train_frac = 0.8,
##                                    min_count_positive_share = 0.3,
##                                    seed = 42L)
##                                    -> list(train, test,
##                                            family_spec,
##                                            scenario_meta,
##                                            strata_summary)
##
## Outcome encoding:
##   y1 = logit(clip(occurrence_seasonal_mean, eps, 1 - eps))
##   y2 = log1p(count_seasonal_mean)
##
## The encoder treats both coordinates as Gaussian on the linear
## predictor scale, matching the homogeneous gaussian family_spec used
## throughout Block 6.

.NE_USA_SUBREGION_MAP <- c(
  "Maine"         = "NE_FOREST",
  "New Hampshire" = "NE_FOREST",
  "Vermont"       = "NE_FOREST",
  "Massachusetts" = "MID_ATL_COAST",
  "Rhode Island"  = "MID_ATL_COAST",
  "Connecticut"   = "MID_ATL_COAST",
  "New Jersey"    = "MID_ATL_COAST",
  "Delaware"      = "MID_ATL_COAST",
  "Maryland"      = "MID_ATL_COAST",
  "Pennsylvania"  = "APPALACHIAN",
  "West Virginia" = "APPALACHIAN",
  "New York"      = "GREAT_LAKES"
)

.assign_state_to_subregion <- function(state_name) {
  unname(.NE_USA_SUBREGION_MAP[state_name])
}

.load_ebird_species_seasonal_p2 <- function(sp, states_sv,
                                              worldclim_rast,
                                              resolution = "9km") {
  occ <- ebirdst::load_raster(species = sp, product = "occurrence",
                                 period = "seasonal",
                                 resolution = resolution)
  cnt <- ebirdst::load_raster(species = sp, product = "count",
                                 period = "seasonal",
                                 resolution = resolution)
  if (!"breeding" %in% names(occ)) {
    stop(sprintf(
      "Species '%s' has no 'breeding' season layer in occurrence raster.",
      sp
    ))
  }
  occ_br <- occ[["breeding"]]
  cnt_br <- cnt[["breeding"]]

  states_ne <- states_sv[states_sv$NAME_1 %in% names(.NE_USA_SUBREGION_MAP), ]
  states_ne_eq <- terra::project(states_ne, terra::crs(occ_br))

  occ_ne <- terra::crop(occ_br, states_ne_eq)
  occ_ne <- terra::mask(occ_ne, states_ne_eq)
  cnt_ne <- terra::crop(cnt_br, states_ne_eq)
  cnt_ne <- terra::mask(cnt_ne, states_ne_eq)

  cell_idx <- which(!is.na(terra::values(occ_ne)) &
                       !is.na(terra::values(cnt_ne)))
  if (length(cell_idx) == 0L) {
    stop(sprintf(
      "No non-NA cells for species '%s' in the NE-USA mask.", sp
    ))
  }
  xy_eq <- terra::xyFromCell(occ_ne, cell_idx)

  cell_pts_eq <- terra::vect(xy_eq, type = "points",
                                crs = terra::crs(occ_ne))
  cell_pts_wgs <- terra::project(cell_pts_eq, "EPSG:4326")
  lonlat <- terra::crds(cell_pts_wgs)

  state_join <- terra::extract(states_ne_eq, cell_pts_eq)
  state_name <- state_join$NAME_1
  sub_region <- .assign_state_to_subregion(state_name)

  occ_vals <- terra::extract(occ_ne, cell_idx)
  cnt_vals <- terra::extract(cnt_ne, cell_idx)
  occ_vec <- if (is.data.frame(occ_vals)) occ_vals[, 1L] else occ_vals
  cnt_vec <- if (is.data.frame(cnt_vals)) cnt_vals[, 1L] else cnt_vals

  bio_vals <- terra::extract(
    worldclim_rast,
    cell_pts_wgs
  )
  bio_cols <- intersect(c("wc2.1_2.5m_bio_1",
                            "wc2.1_2.5m_bio_12",
                            "wc2.1_2.5m_bio_15"),
                         colnames(bio_vals))
  if (length(bio_cols) != 3L) {
    stop(sprintf(
      "WorldClim extract did not return all three layers (got %d): %s",
      length(bio_cols), paste(colnames(bio_vals), collapse = ", ")
    ))
  }
  bio1  <- bio_vals[[bio_cols[1L]]]
  bio12 <- bio_vals[[bio_cols[2L]]]
  bio15 <- bio_vals[[bio_cols[3L]]]

  df <- data.frame(
    species          = sp,
    state            = state_name,
    sub_region       = sub_region,
    lon              = lonlat[, 1L],
    lat              = lonlat[, 2L],
    occurrence_mean  = occ_vec,
    count_mean       = cnt_vec,
    bio1             = bio1,
    bio12            = bio12,
    bio15            = bio15,
    stringsAsFactors = FALSE
  )
  df <- df[stats::complete.cases(df), , drop = FALSE]
  df <- df[!is.na(df$sub_region), , drop = FALSE]
  df
}

.grid_sample_stratified <- function(df, strata_col, n_max, seed = 42L) {
  if (!strata_col %in% colnames(df)) {
    stop(sprintf("Strata column '%s' not found in df.", strata_col))
  }
  set.seed(seed)
  strata_vals <- df[[strata_col]]
  unique_strata <- sort(unique(strata_vals))
  parts <- lapply(unique_strata, function(g) {
    idx <- which(strata_vals == g)
    if (length(idx) <= n_max) return(df[idx, , drop = FALSE])
    keep <- sample.int(length(idx), n_max)
    df[idx[keep], , drop = FALSE]
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

build_ebird_path_a_multi_dataset <- function(
    species = c("reevir1", "amered", "ovenbi1", "btnwar", "woothr"),
    resolution = "9km",
    n_max_per_strata = 750L,
    train_frac = 0.8,
    min_count_positive_share = 0.3,
    eps_logit = 1e-4,
    seed = 42L,
    cache_dir_gadm = file.path("inst", "benchmarks", "cache", "gadm"),
    cache_dir_worldclim = file.path("inst", "benchmarks", "cache",
                                       "worldclim", "climate",
                                       "wc2.1_2.5m")) {
  if (!requireNamespace("ebirdst", quietly = TRUE)) {
    stop("Package 'ebirdst' is required for build_ebird_path_a_multi_dataset().")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for build_ebird_path_a_multi_dataset().")
  }
  if (!requireNamespace("geodata", quietly = TRUE)) {
    stop("Package 'geodata' is required for build_ebird_path_a_multi_dataset().")
  }

  message("Loading GADM USA level-1 polygons ...")
  states_sv <- geodata::gadm(country = "USA", level = 1,
                                path = cache_dir_gadm)

  wc_files <- file.path(cache_dir_worldclim,
                          paste0("wc2.1_2.5m_bio_",
                                  c(1L, 12L, 15L), ".tif"))
  missing_wc <- !file.exists(wc_files)
  if (any(missing_wc)) {
    stop(sprintf(
      "Missing WorldClim layers in cache: %s",
      paste(wc_files[missing_wc], collapse = ", ")
    ))
  }
  message("Stacking WorldClim layers (bio1, bio12, bio15) ...")
  worldclim_rast <- terra::rast(wc_files)

  message(sprintf(
    "Building per-species cell tables (%d species) at resolution %s ...",
    length(species), resolution
  ))
  per_species <- lapply(species, function(sp) {
    message(sprintf("  [%s] loading rasters ...", sp))
    .load_ebird_species_seasonal_p2(sp, states_sv, worldclim_rast,
                                      resolution = resolution)
  })
  full <- do.call(rbind, per_species)
  full$group <- factor(paste(full$species, full$sub_region, sep = "__"))

  strata_summary_raw <- do.call(rbind, lapply(
    split(full, full$group),
    function(d) data.frame(
      group               = as.character(d$group[1L]),
      species             = d$species[1L],
      sub_region          = d$sub_region[1L],
      n_cells             = nrow(d),
      share_count_positive = mean(d$count_mean > 0),
      stringsAsFactors    = FALSE
    )
  ))
  rownames(strata_summary_raw) <- NULL

  keep_groups <- strata_summary_raw$group[
    strata_summary_raw$share_count_positive >=
      min_count_positive_share
  ]
  full_kept <- full[full$group %in% keep_groups, , drop = FALSE]
  full_kept$group <- droplevels(full_kept$group)

  message(sprintf(
    "Keeping %d / %d strata after filter share_count_positive >= %.2f.",
    length(keep_groups), nrow(strata_summary_raw),
    min_count_positive_share
  ))

  sampled <- .grid_sample_stratified(full_kept, strata_col = "group",
                                       n_max = n_max_per_strata,
                                       seed = seed)
  sampled$group <- droplevels(sampled$group)

  sampled$y1 <- stats::qlogis(
    pmin(pmax(sampled$occurrence_mean, eps_logit), 1 - eps_logit)
  )
  sampled$y2 <- log1p(sampled$count_mean)

  set.seed(seed)
  n_total <- nrow(sampled)
  shuffle <- sample.int(n_total)
  sampled <- sampled[shuffle, , drop = FALSE]
  n_train <- floor(train_frac * n_total)
  train_idx <- seq_len(n_train)
  test_idx  <- (n_train + 1L):n_total
  train <- sampled[train_idx, , drop = FALSE]
  test  <- sampled[test_idx,  , drop = FALSE]
  rownames(train) <- NULL
  rownames(test)  <- NULL

  family_spec <- .family_spec_gaussian_homog(p = 2L)

  scenario_meta <- list(
    scenario_id = "ebird_path_a_multi",
    tier        = "external_validation",
    p           = 2L,
    gdpar_amm   = amm_spec(
      p = 2L,
      dims = dimwise(a = ~ lon + lat + bio1 + bio12 + bio15, b = NULL)
    ),
    gdpar_group_formula = ~ group,
    gdpar_rhs   = c("lon", "lat", "bio1", "bio12", "bio15"),
    brms_rhs    =
      "lon + lat + bio1 + bio12 + bio15 + (1 | group)",
    mgcv_formula_per_coord = list(
      stats::as.formula(
        "y1 ~ lon + lat + bio1 + bio12 + bio15 + s(group, bs = \"re\")"
      ),
      stats::as.formula(
        "y2 ~ lon + lat + bio1 + bio12 + bio15 + s(group, bs = \"re\")"
      )
    ),
    inla_formula_per_coord = list(
      stats::as.formula(
        "y1 ~ lon + lat + bio1 + bio12 + bio15 + f(group, model = \"iid\")"
      ),
      stats::as.formula(
        "y2 ~ lon + lat + bio1 + bio12 + bio15 + f(group, model = \"iid\")"
      )
    ),
    rstanarm_formula_per_coord = list(
      stats::as.formula(
        "y1 ~ lon + lat + bio1 + bio12 + bio15 + (1 | group)"
      ),
      stats::as.formula(
        "y2 ~ lon + lat + bio1 + bio12 + bio15 + (1 | group)"
      )
    )
  )

  list(
    train         = train,
    test          = test,
    family_spec   = family_spec,
    scenario_meta = scenario_meta,
    strata_summary = strata_summary_raw,
    keep_groups   = keep_groups
  )
}

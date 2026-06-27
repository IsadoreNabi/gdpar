# Block 9, Axis 2 sub-unit 1: spatial dependence-robust inference (Session
# B9.37, decision D100). Layers 1-3 (weight/Moran algebra, API guards, the
# Moran + permutation statistical property, and the Stan-free RNG-order lock of
# the shared engine) run by default. Layer 4 (end-to-end fit-based behaviour and
# the non-negotiable bit-exact temporal-refactor gate) is gated behind
# GDPAR_RUN_BLOCK9_DEP_FITS=1 because each assertion refits the EB model.

# ---------------------------------------------------------------------------
# Layer 1: spatial weight / Moran's I / resampler algebra (no Stan)
# ---------------------------------------------------------------------------

test_that("spatial default cells-per-axis follows the n^(1/4) variance rate", {
  expect_identical(.gdpar_spatial_default_g(27L), 2L)
  expect_identical(.gdpar_spatial_default_g(100L), 3L)
  expect_identical(.gdpar_spatial_default_g(10000L), 10L)
  expect_true(.gdpar_spatial_default_g(1L) >= 2L)  # floor at 2
})

# ---- data-driven spatial block size (D101) ----

test_that("spatial auto selector returns a valid g and is seed-deterministic", {
  set.seed(1)
  n <- 250L
  co <- cbind(runif(n), runif(n))
  z <- rnorm(n)
  mins <- apply(co, 2L, min)
  ranges <- apply(co, 2L, function(v) diff(range(v)))
  s1 <- .gdpar_spatial_block_length_auto(co, z, "tiled", TRUE, mins, ranges,
                                         seed = 7)
  s2 <- .gdpar_spatial_block_length_auto(co, z, "tiled", TRUE, mins, ranges,
                                         seed = 7)
  expect_true(s1$method %in% c("auto", "rate"))
  expect_gte(s1$block_size, 2L)
  expect_identical(s1$block_size, s2$block_size)   # deterministic given seed
})

test_that("spatial auto selector falls back to the rate on degenerate input", {
  set.seed(2)
  co <- cbind(runif(40L), runif(40L))
  mins <- apply(co, 2L, min)
  ranges <- apply(co, 2L, function(v) diff(range(v)))
  # tiny n -> rate
  small <- .gdpar_spatial_block_length_auto(co[1:18, ], rnorm(18L), "tiled",
                                            TRUE, apply(co[1:18, ], 2L, min),
                                            apply(co[1:18, ], 2L,
                                                  function(v) diff(range(v))),
                                            seed = 7)
  expect_identical(small$method, "rate")
  # zero-variance residual surrogate -> rate
  flat <- .gdpar_spatial_block_length_auto(co, rep(3, 40L), "tiled", TRUE,
                                           mins, ranges, seed = 7)
  expect_identical(flat$method, "rate")
})

test_that("kNN adjacency gives each row exactly k neighbours, no isolation", {
  coords <- as.matrix(expand.grid(x = 1:6, y = 1:6))
  W <- .gdpar_knn_adjacency(coords, 4L)
  expect_true(all(rowSums(W) == 4L))
  expect_true(all(diag(W) == 0))            # no self-neighbour
})

test_that("distance-band adjacency isolates no location", {
  coords <- as.matrix(expand.grid(x = 1:6, y = 1:6))
  W <- .gdpar_distance_band_adjacency(coords)
  expect_true(all(rowSums(W) >= 1L))
})

test_that("Moran's I is negative on a checkerboard and positive on blocks", {
  m <- 8L
  grid <- expand.grid(x = 1:m, y = 1:m)
  coords <- as.matrix(grid)
  Wn <- { W <- .gdpar_knn_adjacency(coords, 4L); W / rowSums(W) }
  checker <- ((grid$x + grid$y) %% 2) * 2 - 1
  blocks  <- ifelse(grid$x <= m / 2, 1, -1)
  expect_lt(.gdpar_morans_i(checker, Wn, S0 = sum(Wn)), 0)
  expect_gt(.gdpar_morans_i(blocks,  Wn, S0 = sum(Wn)), 0)
})

test_that("row-standardized weights satisfy the S0 = n identity", {
  coords <- as.matrix(expand.grid(x = 1:5, y = 1:5))
  Wn <- { W <- .gdpar_knn_adjacency(coords, 3L); W / rowSums(W) }
  expect_equal(sum(Wn), nrow(coords))
})

test_that("spatial block resampler returns exactly n valid indices", {
  coords <- as.matrix(expand.grid(x = 1:8, y = 1:8))
  n <- nrow(coords)
  mins <- apply(coords, 2L, min)
  ranges <- apply(coords, 2L, function(z) diff(range(z)))
  for (scheme in c("tiled", "moving")) {
    set.seed(3)
    idx <- .gdpar_spatial_block_indices(coords, 3L, scheme, TRUE, mins, ranges)
    expect_length(idx, n)
    expect_true(is.integer(idx))
    expect_true(all(idx >= 1L & idx <= n))
  }
})

test_that("spatial block resampler is deterministic under a fixed seed", {
  coords <- as.matrix(expand.grid(x = 1:8, y = 1:8))
  mins <- apply(coords, 2L, min)
  ranges <- apply(coords, 2L, function(z) diff(range(z)))
  set.seed(9)
  a <- .gdpar_spatial_block_indices(coords, 3L, "tiled", TRUE, mins, ranges)
  set.seed(9)
  b <- .gdpar_spatial_block_indices(coords, 3L, "tiled", TRUE, mins, ranges)
  expect_identical(a, b)
})

test_that("GOLDEN: shared engine RNG-order is frozen (Stan-free, D100)", {
  # The engine draws the B per-refit Stan seeds FIRST, then calls resample_fun()
  # once per iteration. Freezing the post-seed-draw resample sequence locks that
  # order: reordering it (resample before the seed draw) trips this golden,
  # which is the Stan-free half of the bit-exact refactor gate.
  n <- 30L; bl <- 3L; B <- 3L; type <- "moving"
  set.seed(2026L)
  refit_seeds <- sample.int(.Machine$integer.max, B)
  first  <- .gdpar_block_bootstrap_data_indices(n, bl, type)
  second <- .gdpar_block_bootstrap_data_indices(n, bl, type)
  expect_identical(refit_seeds, c(853315193L, 601901350L, 237766255L))
  expect_identical(
    first,
    c(27L, 28L, 29L, 12L, 13L, 14L, 4L, 5L, 6L, 16L, 17L, 18L, 5L, 6L,
      7L, 12L, 13L, 14L, 26L, 27L, 28L, 2L, 3L, 4L, 24L, 25L, 26L, 19L,
      20L, 21L)
  )
  expect_identical(
    second,
    c(15L, 16L, 17L, 14L, 15L, 16L, 17L, 18L, 19L, 12L, 13L, 14L, 10L,
      11L, 12L, 22L, 23L, 24L, 18L, 19L, 20L, 2L, 3L, 4L, 5L, 6L, 7L,
      12L, 13L, 14L)
  )
})

# ---------------------------------------------------------------------------
# Layer 2: public API guards (no Stan; fail before any refit)
# ---------------------------------------------------------------------------

fake_eb <- structure(
  list(conditional_fit = list(), family = list(name = "gaussian")),
  class = c("gdpar_eb_fit", "list")
)

test_that("spatial functions validate the supplied fit (EB or FB, D102)", {
  # A non-fit object is rejected outright.
  expect_error(
    gdpar_spatial_dependence_diagnostic(list(1, 2, 3),
                                        coords = cbind(1:3, 1:3)),
    class = "gdpar_input_error"
  )
  # A scalar full-Bayes gdpar_fit is accepted as a class (D102), but a malformed
  # one lacking its HMC fit is caught.
  no_fit <- structure(list(), class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_spatial_dependence_diagnostic(no_fit, coords = cbind(1:3, 1:3)),
    class = "gdpar_internal_error"
  )
  # Multivariate full-Bayes fits remain deferred.
  multi_fb <- structure(list(fit = list(), p = 2L),
                        class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_spatial_dependence_robust(multi_fb, data = data.frame(x = 1),
                                    coords = cbind(1, 1)),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("spatial robust rejects a non-data-frame data argument", {
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = matrix(1:4, 2),
                                    coords = cbind(1:2, 1:2)),
    class = "gdpar_input_error"
  )
})

test_that("coordinate validation rejects malformed coords", {
  df <- data.frame(x = rnorm(10), y = rnorm(10))
  # wrong number of columns
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = df,
                                    coords = matrix(rnorm(30), 10, 3)),
    class = "gdpar_input_error"
  )
  # wrong number of rows
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = df,
                                    coords = matrix(rnorm(10), 5, 2)),
    class = "gdpar_input_error"
  )
  # non-finite coordinate
  bad <- cbind(c(rnorm(9), NA), rnorm(10))
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = df, coords = bad),
    class = "gdpar_input_error"
  )
})

test_that("spatial robust aborts on collinear coordinates", {
  df <- data.frame(x = rnorm(10), y = rnorm(10))
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = df,
                                    coords = cbind(1:10, rep(5, 10))),
    class = "gdpar_input_error"
  )
})

test_that("spatial robust rejects a non-positive block_size", {
  df <- data.frame(x = rnorm(10), y = rnorm(10))
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = df,
                                    coords = cbind(runif(10), runif(10)),
                                    block_size = 0L),
    class = "gdpar_input_error"
  )
})

test_that("spatial robust rejects a non-\"auto\" string block_size", {
  # Caught Stan-free, before coords validation / require_suggested / refit.
  expect_error(
    gdpar_spatial_dependence_robust(fake_eb, data = data.frame(x = 1, y = 1),
                                    coords = cbind(1, 1), block_size = "bogus"),
    class = "gdpar_input_error"
  )
})

# ---------------------------------------------------------------------------
# Layer 3: statistical detection property (Moran + permutation; no Stan)
# ---------------------------------------------------------------------------

test_that("Moran + permutation flags spatial clustering and not iid noise", {
  set.seed(11)
  m <- 12L
  grid <- expand.grid(x = 1:m, y = 1:m)
  coords <- as.matrix(grid)
  n <- nrow(coords)
  Wn <- { W <- .gdpar_knn_adjacency(coords, 6L); W / rowSums(W) }
  S0 <- sum(Wn)
  EI <- -1 / (n - 1)
  perm_p <- function(field) {
    I <- .gdpar_morans_i(field, Wn, S0 = S0)
    dev <- abs(I - EI)
    ge <- sum(replicate(
      199L, abs(.gdpar_morans_i(sample(field), Wn, S0 = S0) - EI) >= dev
    ))
    list(I = I, p = (1 + ge) / 200)
  }
  clustered <- as.numeric(grid$x)        # smooth gradient: strong positive I
  noise <- rnorm(n)
  rc <- perm_p(clustered)
  rn <- perm_p(noise)
  expect_gt(rc$I, 0.3)
  expect_lt(rc$p, 0.05)
  expect_gt(rn$p, 0.05)
})

test_that("spatial auto adapts: shorter dependence range picks smaller blocks", {
  # A Gaussian random field with exponential correlation of range rho. Shorter
  # range (weaker, more local dependence) should be captured by smaller blocks
  # (larger g) than a long-range, strongly dependent field. The corrected
  # bias anchor + 1/n_tiles variance term make the calibration adaptive (D101);
  # the inverted-anchor / jackknife-MC form would collapse onto g_min always.
  field_g <- function(rho) {
    set.seed(21L)                      # same coords + noise; only rho differs
    co <- cbind(runif(300L), runif(300L))
    L <- chol(exp(-as.matrix(stats::dist(co)) / rho) + diag(1e-6, 300L))
    z <- as.numeric(t(L) %*% rnorm(300L))
    mins <- apply(co, 2L, min)
    ranges <- apply(co, 2L, function(v) diff(range(v)))
    .gdpar_spatial_block_length_auto(co, z, "tiled", TRUE, mins, ranges,
                                     seed = 7)
  }
  g_long <- field_g(0.40)
  g_short <- field_g(0.05)
  expect_identical(g_long$method, "auto")
  expect_identical(g_short$method, "auto")
  expect_gte(g_short$block_size, g_long$block_size)
})

# ---------------------------------------------------------------------------
# Layer 4: end-to-end behaviour + bit-exact refactor gate (gated; refits)
# ---------------------------------------------------------------------------

test_that("spatial diagnostic detects a spatial trend and not iid residuals", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(540)
  n <- 144L
  gx <- runif(n); gy <- runif(n)
  x <- rnorm(n)
  y_sp  <- 1 + 0.5 * x + 3 * (gx + gy) + rnorm(n, sd = 0.3)  # omitted spatial trend
  y_iid <- 1 + 0.5 * x + rnorm(n)
  df_sp  <- data.frame(x = x, y = y_sp)
  df_iid <- data.frame(x = x, y = y_iid)
  fit_sp <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df_sp,
                     chains = 2, iter_warmup = 150, iter_sampling = 150,
                     verbose = FALSE, refresh = 0, seed = 54)
  fit_iid <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df_iid,
                      chains = 2, iter_warmup = 150, iter_sampling = 150,
                      verbose = FALSE, refresh = 0, seed = 55)
  coords <- cbind(gx, gy)
  d_sp <- gdpar_spatial_dependence_diagnostic(fit_sp, coords = coords,
                                              residual_type = "response",
                                              seed = 1)
  d_iid <- gdpar_spatial_dependence_diagnostic(fit_iid, coords = coords,
                                               residual_type = "response",
                                               seed = 1)
  expect_s3_class(d_sp, "gdpar_spatial_dependence_diagnostic")
  expect_gt(d_sp$morans_i, 0.1)
  expect_lt(d_sp$p_value, 0.05)
  expect_gt(d_iid$p_value, 0.05)
})

test_that("spatial robust returns a robust-SE table and is seed-deterministic", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(560)
  n <- 120L
  gx <- runif(n); gy <- runif(n)
  x <- rnorm(n)
  y <- 1 + 0.5 * x + 3 * (gx + gy) + rnorm(n, sd = 0.3)
  df <- data.frame(x = x, y = y)
  fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
                  chains = 2, iter_warmup = 150, iter_sampling = 150,
                  verbose = FALSE, refresh = 0, seed = 56)
  coords <- cbind(gx, gy)
  rb <- gdpar_spatial_dependence_robust(fit, data = df, coords = coords,
                                        B = 12, seed = 7, verbose = FALSE,
                                        iter_warmup = 150, iter_sampling = 150,
                                        chains = 2)
  expect_s3_class(rb, "gdpar_spatial_dependence_robust")
  expect_true(all(c("estimate", "model_se", "robust_se", "se_ratio",
                    "ci_lower", "ci_upper") %in% names(rb$table)))
  expect_true(all(is.finite(rb$table$robust_se)))
  expect_gte(rb$B_ok, 2L)
  expect_true(all(rb$table$ci_lower <= rb$table$ci_upper))
  expect_gt(max(rb$table$se_ratio), 1)            # dependence inflates the SE
  rb2 <- gdpar_spatial_dependence_robust(fit, data = df, coords = coords,
                                         B = 12, seed = 7, verbose = FALSE,
                                         iter_warmup = 150, iter_sampling = 150,
                                         chains = 2)
  expect_equal(rb$table$robust_se, rb2$table$robust_se)
})

test_that("block_size = \"auto\" selects a data-driven g end-to-end (D101)", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(561)
  n <- 150L
  gx <- runif(n); gy <- runif(n)
  x <- rnorm(n)
  y <- 1 + 0.5 * x + 3 * (gx + gy) + rnorm(n, sd = 0.3)
  df <- data.frame(x = x, y = y)
  fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
                  chains = 2, iter_warmup = 150, iter_sampling = 150,
                  verbose = FALSE, refresh = 0, seed = 57)
  coords <- cbind(gx, gy)
  rb <- gdpar_spatial_dependence_robust(fit, data = df, coords = coords,
                                        block_size = "auto", B = 12, seed = 7,
                                        verbose = FALSE, iter_warmup = 150,
                                        iter_sampling = 150, chains = 2)
  expect_s3_class(rb, "gdpar_spatial_dependence_robust")
  expect_true(rb$block_size_method %in% c("auto", "rate"))
  expect_gte(rb$block_size, 2L)
  expect_true(all(is.finite(rb$table$robust_se)))
  # The auto g is a deterministic function of the (seeded) residual field.
  rb2 <- gdpar_spatial_dependence_robust(fit, data = df, coords = coords,
                                         block_size = "auto", B = 12, seed = 7,
                                         verbose = FALSE, iter_warmup = 150,
                                         iter_sampling = 150, chains = 2)
  expect_identical(rb$block_size, rb2$block_size)
})

test_that("REFACTOR GATE: temporal path is bit-identical pre/post engine (D100)", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(770)
  n <- 90L
  x <- rnorm(n)
  y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.5), n))
  df <- data.frame(x = x, y = y, t = seq_len(n))
  fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
                  chains = 2, iter_warmup = 150, iter_sampling = 150,
                  verbose = FALSE, refresh = 0, seed = 77)
  args <- list(object = fit, data = df, index = df$t, block_length = 4L,
               type = "moving", B = 8, level = 0.95, seed = 123,
               verbose = FALSE, iter_warmup = 150, iter_sampling = 150,
               chains = 2)
  new_res <- do.call(gdpar_dependence_robust, args)
  old_res <- do.call(.gdpar_dependence_robust_pre_refactor, args)
  expect_identical(new_res$table, old_res$table)
  expect_identical(new_res$block_length, old_res$block_length)
  expect_identical(new_res$B_ok, old_res$B_ok)
  expect_identical(new_res$seed, old_res$seed)
})

# ---------------------------------------------------------------------------
# Layer 3b: full-Bayes (gdpar_fit) spatial end-to-end (gated; D102)
# ---------------------------------------------------------------------------

test_that("FB path: spatial diagnostic + robust SE end-to-end (D102)", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(580)
  n <- 120L
  gx <- runif(n); gy <- runif(n)
  x <- rnorm(n)
  y <- 1 + 0.5 * x + 3 * (gx + gy) + rnorm(n, sd = 0.3)
  df <- data.frame(x = x, y = y)
  coords <- cbind(gx, gy)
  fb <- gdpar(y ~ x, amm = amm_spec(a = ~ x), data = df,
              chains = 2, iter_warmup = 150, iter_sampling = 150,
              verbose = FALSE, refresh = 0, seed = 58)
  expect_s3_class(fb, "gdpar_fit")
  # The same spatial diagnostic accepts a gdpar_fit (D102).
  d <- gdpar_spatial_dependence_diagnostic(fb, coords = coords, seed = 1)
  expect_s3_class(d, "gdpar_spatial_dependence_diagnostic")
  # Spatial robust SE by full-Bayes refit. B kept tiny for the gated suite.
  rb <- gdpar_spatial_dependence_robust(fb, data = df, coords = coords,
                                        B = 4, seed = 7, verbose = FALSE,
                                        iter_warmup = 150, iter_sampling = 150,
                                        chains = 2)
  expect_s3_class(rb, "gdpar_spatial_dependence_robust")
  expect_true(all(c("estimate", "model_se", "robust_se", "se_ratio",
                    "ci_lower", "ci_upper") %in% names(rb$table)))
  expect_true(all(is.finite(rb$table$robust_se)))
  expect_gte(rb$B_ok, 2L)
  expect_false(is.null(rb$refit_diagnostics))
})

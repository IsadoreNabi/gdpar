# Sub-phase 8.3.9 (2026-05-22): structural tests for the unified
# bootstrap helpers + manifest CSV scaffolding (R/golden_helpers.R).
#
# The tests do not bootstrap any actual Stan fits; they exercise the
# helper functions (path resolution, init / read / add idempotency,
# canonical column set, roster, fit_code_hash, version detection,
# RDS filename convention). Bootstrap of the 14 goldens is gated by
# per-config env vars in test-golden_regression_K2.R and similar.

test_that(".gdpar_golden_manifest_columns returns the canonical 16 columns", {
  cols <- .gdpar_golden_manifest_columns()
  expect_length(cols, 16L)
  for (key in c("family", "K", "basis_type", "p", "sub_phase",
                "bootstrap_date", "cmdstan_version", "R_version",
                "DHARMa_version", "fit_code_hash", "seed",
                "n_chains", "n_iter_warmup", "n_iter_sampling",
                "rds_path", "n_obs")) {
    expect_true(key %in% cols, info = sprintf("missing column %s", key))
  }
})

test_that(".gdpar_golden_manifest_path returns string", {
  p <- .gdpar_golden_manifest_path()
  expect_type(p, "character")
  expect_match(p, "golden_manifest.csv$")
})

test_that(".gdpar_golden_manifest_init creates empty CSV if missing", {
  tmp <- file.path(tempdir(), paste0("gm_", format(Sys.time(), "%H%M%S%OS3"),
                                     ".csv"))
  on.exit(unlink(tmp))
  .gdpar_golden_manifest_init(tmp)
  expect_true(file.exists(tmp))
  df <- .gdpar_golden_manifest_read(tmp)
  expect_equal(nrow(df), 0L)
  expect_true(all(.gdpar_golden_manifest_columns() %in% names(df)))
})

test_that(".gdpar_golden_manifest_init is idempotent (no-op if exists)", {
  tmp <- file.path(tempdir(), paste0("gm_", format(Sys.time(), "%H%M%S%OS3"),
                                     ".csv"))
  on.exit(unlink(tmp))
  .gdpar_golden_manifest_init(tmp)
  before <- readLines(tmp)
  .gdpar_golden_manifest_init(tmp)
  expect_identical(readLines(tmp), before)
})

test_that(".gdpar_golden_manifest_add appends a row", {
  tmp <- file.path(tempdir(), paste0("gm_", format(Sys.time(), "%H%M%S%OS3"),
                                     ".csv"))
  on.exit(unlink(tmp))
  .gdpar_golden_manifest_add(
    list(family = "tweedie", K = 3L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.5b", bootstrap_date = "2026-05-22",
         cmdstan_version = "2.38.0", R_version = "4.6.0",
         fit_code_hash = "deadbeef", seed = 42L, n_chains = 1L,
         n_iter_warmup = 200L, n_iter_sampling = 200L,
         rds_path = "data/golden_K3_tweedie_K3.rds", n_obs = 80L),
    path = tmp
  )
  df <- .gdpar_golden_manifest_read(tmp)
  expect_equal(nrow(df), 1L)
  expect_equal(df$family[1L], "tweedie")
  expect_equal(as.character(df$K[1L]), "3")
})

test_that(".gdpar_golden_manifest_add is idempotent on (family,K,basis_type,p)", {
  tmp <- file.path(tempdir(), paste0("gm_", format(Sys.time(), "%H%M%S%OS3"),
                                     ".csv"))
  on.exit(unlink(tmp))
  base <- list(family = "tweedie", K = 3L, basis_type = "polynomial",
               p = 1L, bootstrap_date = "2026-05-22", seed = 42L)
  .gdpar_golden_manifest_add(base, path = tmp)
  base$bootstrap_date <- "2026-05-23"
  base$seed <- 43L
  .gdpar_golden_manifest_add(base, path = tmp)
  df <- .gdpar_golden_manifest_read(tmp)
  expect_equal(nrow(df), 1L)
  expect_equal(df$bootstrap_date[1L], "2026-05-23")
  expect_equal(as.character(df$seed[1L]), "43")
})

test_that(".gdpar_golden_manifest_add rejects rows lacking required keys", {
  tmp <- file.path(tempdir(), paste0("gm_", format(Sys.time(), "%H%M%S%OS3"),
                                     ".csv"))
  on.exit(unlink(tmp))
  expect_error(
    .gdpar_golden_manifest_add(list(family = "tweedie"), path = tmp),
    class = "gdpar_internal_error"
  )
})

test_that(".gdpar_golden_rds_name follows convention", {
  expect_equal(
    .gdpar_golden_rds_name("beta", K = 2L),
    "golden_K2_beta_K2.rds"
  )
  expect_equal(
    .gdpar_golden_rds_name("gaussian", K = 1L, basis_type = "bspline"),
    "golden_K1_gaussian_K1_bspline.rds"
  )
  expect_equal(
    .gdpar_golden_rds_name("gaussian", K = 1L, basis_type = "polynomial",
                           p = 2L),
    "golden_K1_gaussian_K1_p2.rds"
  )
  expect_equal(
    .gdpar_golden_rds_name("gaussian", K = 1L, basis_type = "bspline",
                           p = 2L),
    "golden_K1_gaussian_K1_bspline_p2.rds"
  )
})

test_that(".gdpar_golden_roster_8_3_9 returns the 14-config roster (Decision D1=1D)", {
  r <- .gdpar_golden_roster_8_3_9()
  expect_length(r, 14L)
  fams <- vapply(r, function(x) x$family, character(1L))
  Ks   <- vapply(r, function(x) as.integer(x$K), integer(1L))
  bts  <- vapply(r, function(x) x$basis_type, character(1L))
  expect_true("tweedie" %in% fams)
  expect_true("student_t" %in% fams)
  expect_true("zip" %in% fams)
  expect_true("zinb" %in% fams)
  expect_true("hurdle_poisson" %in% fams)
  expect_true("hurdle_neg_binomial_2" %in% fams)
  expect_true("lognormal_loc_scale" %in% fams)
  expect_true("het_gauss_beta" %in% fams)
  expect_true("het_gauss_gamma" %in% fams)
  expect_true("het_nb_beta" %in% fams)
  expect_true("bspline" %in% bts)
})

test_that(".gdpar_golden_fit_code_hash is deterministic on a closure", {
  f <- function() {
    set.seed(42L)
    rnorm(10L)
  }
  h1 <- .gdpar_golden_fit_code_hash(f)
  h2 <- .gdpar_golden_fit_code_hash(f)
  expect_identical(h1, h2)
  expect_type(h1, "character")
})

test_that(".gdpar_golden_fit_code_hash distinguishes different closures", {
  f1 <- function() rnorm(10L)
  f2 <- function() rnorm(20L)
  h1 <- .gdpar_golden_fit_code_hash(f1)
  h2 <- .gdpar_golden_fit_code_hash(f2)
  expect_false(identical(h1, h2))
})

test_that(".gdpar_cmdstan_version returns string or NA", {
  v <- .gdpar_cmdstan_version()
  expect_true(is.character(v))
  expect_length(v, 1L)
})

test_that(".gdpar_dharma_version returns string or NA", {
  v <- .gdpar_dharma_version()
  expect_true(is.character(v))
  expect_length(v, 1L)
})

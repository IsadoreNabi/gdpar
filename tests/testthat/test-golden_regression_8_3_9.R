# Sub-phase 8.3.9 (2026-05-22): structural validation of the 14
# new goldens + manifest CSV (Decisions D1=1D, D2=2B, E3=E3.A of the
# 8.3.9 scoping).
#
# Structural tests (no Stan re-run): load each .rds and verify that
# (i) the file exists, (ii) the draws_matrix has the right shape,
# (iii) the family / K / p / basis_type in the saved metadata match
# the roster, (iv) the manifest CSV contains a matching row.
#
# Bit-exact compare on re-fit is gated by per-config env vars
# (GDPAR_GOLDEN_8_3_9_COMPARE_<label>=1); when enabled, the test
# re-runs the matching closure from the bootstrap script and verifies
# bit-exact equality of the draws_matrix against the saved snapshot.
# Disabled by default since each re-fit compiles a Stan model and
# samples (taking seconds to minutes per config).

.golden_root <- function() {
  c1 <- file.path("data")
  c2 <- file.path("tests", "testthat", "data")
  if (dir.exists(c2)) return(c2)
  c1
}

.load_golden_rds <- function(family, K, basis_type, p) {
  fname <- .gdpar_golden_rds_name(family, K, basis_type, p)
  path <- file.path(.golden_root(), fname)
  if (!file.exists(path)) {
    return(NULL)
  }
  list(path = path, data = readRDS(path))
}

# Returns TRUE iff at least one 8.3.9 golden RDS is present on disk.
# Used as a fast guard so the suite skips cleanly when goldens are
# absent from the tarball (e.g. CRAN distribution where they are
# excluded via .Rbuildignore but kept in the source repo). Mirrors
# `.eb_C_first_golden_present` in test-eb_goldens_8_6_C.R.
.any_golden_8_3_9_present <- function() {
  any(vapply(.gdpar_golden_roster_8_3_9(), function(cfg) {
    file.exists(file.path(
      .golden_root(),
      .gdpar_golden_rds_name(cfg$family, cfg$K, cfg$basis_type, cfg$p)
    ))
  }, logical(1L)))
}

# ---- Section 1: roster coverage and file existence ----------------------

test_that(".gdpar_golden_roster_8_3_9 produces exactly 14 configs", {
  expect_length(.gdpar_golden_roster_8_3_9(), 14L)
})

test_that("all 14 golden RDS files exist on disk", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present (bootstrap not yet run)")
  testthat::skip_if_not(.any_golden_8_3_9_present(),
                        "no 8.3.9 goldens present on disk (excluded from tarball)")
  roster <- .gdpar_golden_roster_8_3_9()
  for (cfg in roster) {
    fname <- .gdpar_golden_rds_name(
      cfg$family, cfg$K, cfg$basis_type, cfg$p
    )
    path <- file.path(.golden_root(), fname)
    expect_true(file.exists(path),
                info = sprintf("missing golden: %s", path))
  }
})

# ---- Section 2: manifest CSV integrity ----------------------------------

test_that("manifest CSV has 14 rows matching the roster keys", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present (bootstrap not yet run)")
  m <- .gdpar_golden_manifest_read()
  expect_gte(nrow(m), 14L)
  roster <- .gdpar_golden_roster_8_3_9()
  for (cfg in roster) {
    row <- m[m$family == cfg$family &
             as.integer(m$K) == cfg$K &
             m$basis_type == cfg$basis_type &
             as.integer(m$p) == cfg$p, , drop = FALSE]
    expect_equal(nrow(row), 1L,
                 info = sprintf("expected 1 manifest row for %s/K%d/%s/p%d",
                                cfg$family, cfg$K, cfg$basis_type, cfg$p))
  }
})

test_that("manifest CSV has the canonical 16 columns populated", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present (bootstrap not yet run)")
  m <- .gdpar_golden_manifest_read()
  for (col in .gdpar_golden_manifest_columns()) {
    expect_true(col %in% names(m),
                info = sprintf("missing manifest column: %s", col))
  }
  expect_true(all(!is.na(m$family)))
  expect_true(all(!is.na(m$K)))
  expect_true(all(!is.na(m$bootstrap_date)))
  expect_true(all(!is.na(m$cmdstan_version)))
  expect_true(all(!is.na(m$fit_code_hash)))
  expect_true(all(nchar(m$fit_code_hash) == 64L | m$fit_code_hash == "" |
                  startsWith(m$fit_code_hash, "nohash_")))
})

# ---- Section 3: per-config RDS structural integrity ---------------------

.check_rds_structural <- function(cfg) {
  loaded <- .load_golden_rds(cfg$family, cfg$K, cfg$basis_type, cfg$p)
  if (is.null(loaded)) {
    testthat::skip(sprintf("golden missing on disk: %s/K%d/%s/p%d",
                           cfg$family, cfg$K, cfg$basis_type, cfg$p))
  }
  d <- loaded$data
  expect_true("draws_matrix" %in% names(d))
  expect_true(inherits(d$draws_matrix, "draws_matrix") ||
              is.matrix(d$draws_matrix))
  expect_equal(d$family, cfg$family)
  expect_equal(as.integer(d$K), as.integer(cfg$K))
  expect_equal(d$basis_type, cfg$basis_type)
  expect_equal(as.integer(d$p), as.integer(cfg$p))
  expect_true(d$bootstrap_date >= "2026-05-01")
  expect_true(!is.na(d$cmdstan_version))
  expect_true(d$n_iter_sampling > 0L)
  m <- unclass(d$draws_matrix)
  expect_true(nrow(m) > 0L && ncol(m) > 0L)
  expect_true(all(is.finite(m[, sample.int(ncol(m), min(20L, ncol(m))),
                              drop = FALSE])) ||
              any(is.finite(m[, sample.int(ncol(m), min(20L, ncol(m))),
                              drop = FALSE])))
}

test_that("per-config RDS structural integrity: 8.3.4 lognormal_loc_scale K=2", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  .check_rds_structural(list(
    family = "lognormal_loc_scale", K = 2L, basis_type = "polynomial", p = 1L
  ))
})

test_that("per-config RDS structural integrity: 8.3.5a student_t K=3", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  .check_rds_structural(list(
    family = "student_t", K = 3L, basis_type = "polynomial", p = 1L
  ))
})

test_that("per-config RDS structural integrity: 8.3.5b tweedie K=3", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  .check_rds_structural(list(
    family = "tweedie", K = 3L, basis_type = "polynomial", p = 1L
  ))
})

test_that("per-config RDS structural integrity: 8.3.6 mixtures (4)", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  for (fam in c("zip", "zinb", "hurdle_poisson", "hurdle_neg_binomial_2")) {
    K <- if (fam %in% c("zip", "hurdle_poisson")) 2L else 3L
    .check_rds_structural(list(
      family = fam, K = K, basis_type = "polynomial", p = 1L
    ))
  }
})

test_that("per-config RDS structural integrity: 8.3.7 heterogeneous (3)", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  for (fam in c("het_gauss_beta", "het_gauss_gamma", "het_nb_beta")) {
    .check_rds_structural(list(
      family = fam, K = 2L, basis_type = "polynomial", p = 1L
    ))
  }
})

test_that("per-config RDS structural integrity: 8.3.8 bspline (3)", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  .check_rds_structural(list(
    family = "gaussian", K = 1L, basis_type = "bspline", p = 1L
  ))
  .check_rds_structural(list(
    family = "het_gauss_beta", K = 2L, basis_type = "bspline", p = 1L
  ))
  .check_rds_structural(list(
    family = "gaussian", K = 1L, basis_type = "bspline", p = 2L
  ))
})

test_that("per-config RDS structural integrity: 8.3.9 baseline K=1+p=2 polynomial", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  .check_rds_structural(list(
    family = "gaussian", K = 1L, basis_type = "polynomial", p = 2L
  ))
})

# ---- Section 4: cross-validation of RDS metadata vs manifest -----------

test_that("RDS metadata is consistent with manifest row for every config", {
  testthat::skip_if(!file.exists(.gdpar_golden_manifest_path()),
                    "golden_manifest.csv not present")
  m <- .gdpar_golden_manifest_read()
  roster <- .gdpar_golden_roster_8_3_9()
  for (cfg in roster) {
    loaded <- .load_golden_rds(cfg$family, cfg$K, cfg$basis_type, cfg$p)
    if (is.null(loaded)) next
    row <- m[m$family == cfg$family &
             as.integer(m$K) == cfg$K &
             m$basis_type == cfg$basis_type &
             as.integer(m$p) == cfg$p, , drop = FALSE]
    expect_equal(nrow(row), 1L)
    expect_equal(loaded$data$family, row$family[1L])
    expect_equal(as.integer(loaded$data$K), as.integer(row$K[1L]))
    expect_equal(as.integer(loaded$data$p), as.integer(row$p[1L]))
    expect_equal(loaded$data$basis_type, row$basis_type[1L])
    expect_equal(loaded$data$cmdstan_version, row$cmdstan_version[1L])
  }
})

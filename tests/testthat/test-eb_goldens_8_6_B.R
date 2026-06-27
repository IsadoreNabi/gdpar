# Sub-phase 8.6.B (2026-05-25): structural validation of the 17 EB
# goldens + manifest CSV rows (Charter §3.2 / parcial §5.2 of
# HANDOFF_SUBFASE_8_6_B_PARCIAL.md). Mirror of
# test-golden_regression_8_3_9.R for the EB sub-phase.
#
# Structural tests (no Stan re-run): load each .rds and verify that
# (i) the file exists, (ii) the recorded family / regime / status
# match the roster, (iii) the manifest CSV contains a matching row
# with sub_phase = "8.6.B", (iv) the recorded numerics are finite for
# fitable cases and the guard_class is one of the canonical error
# classes for guard cases. Bootstrap of the 17 goldens lives in
# data-raw/bootstrap_eb_goldens.R and is gated by
# GDPAR_BOOTSTRAP_EB_GOLDENS=1.

.eb_golden_root <- function() {
  c2 <- file.path("tests", "testthat", "data")
  if (dir.exists(c2)) return(c2)
  file.path("data")
}

.load_eb_golden_rds <- function(family, regime) {
  fname <- .gdpar_golden_eb_rds_name(family, regime)
  path <- file.path(.eb_golden_root(), fname)
  if (!file.exists(path)) return(NULL)
  list(path = path, data = readRDS(path))
}

# ---- Section 1: roster coverage and naming -------------------------------

test_that(".gdpar_golden_roster_8_6_B returns exactly 17 entries", {
  expect_length(.gdpar_golden_roster_8_6_B(), 17L)
})

test_that(".gdpar_golden_roster_8_6_B has 5 fitable + 12 guard cases", {
  roster <- .gdpar_golden_roster_8_6_B()
  statuses <- vapply(roster, function(x) x$status, character(1L))
  expect_equal(sum(statuses == "fitable"), 5L)
  expect_equal(sum(statuses == "guard"), 12L)
})

test_that(".gdpar_golden_eb_rds_name follows the canonical convention", {
  expect_equal(
    .gdpar_golden_eb_rds_name("gaussian", "poly"),
    "golden_eb_gaussian_poly.rds"
  )
  expect_equal(
    .gdpar_golden_eb_rds_name("gaussian", "bspline"),
    "golden_eb_gaussian_bspline.rds"
  )
  expect_equal(
    .gdpar_golden_eb_rds_name("zip", "guard"),
    "golden_eb_zip_guard.rds"
  )
})

test_that("EB roster guard_class is canonical for every guard entry", {
  roster <- .gdpar_golden_roster_8_6_B()
  allowed <- c("gdpar_unsupported_feature_error",
               "gdpar_input_error")
  for (cfg in roster) {
    if (cfg$status == "guard") {
      expect_true(cfg$guard_class %in% allowed,
                  info = sprintf("guard_class %s not canonical for %s/%s",
                                 cfg$guard_class, cfg$family, cfg$regime))
    } else {
      expect_true(is.na(cfg$guard_class),
                  info = sprintf("fitable %s/%s should have NA guard_class",
                                 cfg$family, cfg$regime))
    }
  }
})

# ---- Section 2: golden RDS file presence (skip if not bootstrapped) -----

.eb_first_golden_present <- function() {
  any(vapply(.gdpar_golden_roster_8_6_B(), function(cfg) {
    file.exists(file.path(.eb_golden_root(),
                          .gdpar_golden_eb_rds_name(cfg$family, cfg$regime)))
  }, logical(1L)))
}

test_that("all 17 EB golden RDS files exist on disk (skip if not bootstrapped)", {
  testthat::skip_if_not(.eb_first_golden_present(),
                        "no EB goldens present (bootstrap not yet run)")
  for (cfg in .gdpar_golden_roster_8_6_B()) {
    path <- file.path(.eb_golden_root(),
                      .gdpar_golden_eb_rds_name(cfg$family, cfg$regime))
    expect_true(file.exists(path),
                info = sprintf("missing golden: %s", path))
  }
})

# ---- Section 3: per-config RDS structural integrity (when present) ------

test_that("fitable EB goldens carry finite theta_ref_hat / theta_ref_se / kappa", {
  testthat::skip_if_not(.eb_first_golden_present(),
                        "no EB goldens present")
  for (cfg in .gdpar_golden_roster_8_6_B()) {
    if (cfg$status != "fitable") next
    loaded <- .load_eb_golden_rds(cfg$family, cfg$regime)
    if (is.null(loaded)) {
      testthat::skip(sprintf("not bootstrapped: %s/%s",
                             cfg$family, cfg$regime))
    }
    d <- loaded$data
    expect_equal(d$family, cfg$family)
    expect_equal(d$regime, cfg$regime)
    expect_equal(d$status, "fitable")
    expect_equal(d$sub_phase, "8.6.B")
    expect_true(all(is.finite(d$theta_ref_hat)),
                info = sprintf("non-finite theta_ref_hat in %s/%s",
                               cfg$family, cfg$regime))
    expect_true(all(is.finite(d$theta_ref_se)),
                info = sprintf("non-finite theta_ref_se in %s/%s",
                               cfg$family, cfg$regime))
    expect_true(is.finite(d$diagnostics_numerical$kappa))
    expect_true(d$correction_applied)
  }
})

test_that("guard EB goldens carry canonical guard_class metadata", {
  testthat::skip_if_not(.eb_first_golden_present(),
                        "no EB goldens present")
  for (cfg in .gdpar_golden_roster_8_6_B()) {
    if (cfg$status != "guard") next
    loaded <- .load_eb_golden_rds(cfg$family, cfg$regime)
    if (is.null(loaded)) next
    d <- loaded$data
    expect_equal(d$status, "guard")
    expect_equal(d$family, cfg$family)
    expect_equal(d$guard_class, cfg$guard_class)
  }
})

# ---- Section 4: manifest CSV cross-validation ----------------------------

test_that("manifest CSV contains sub_phase = '8.6.B' rows when bootstrapped", {
  testthat::skip_if_not(file.exists(.gdpar_golden_manifest_path()),
                        "golden_manifest.csv not present")
  testthat::skip_if_not(.eb_first_golden_present(),
                        "no EB goldens present")
  m <- .gdpar_golden_manifest_read()
  eb_rows <- m[m$sub_phase == "8.6.B", , drop = FALSE]
  expect_true(nrow(eb_rows) >= 1L)
})

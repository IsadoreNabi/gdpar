# Sub-phase 8.6.C (Session 12 closure, 2026-05-25): structural
# validation of the 17 EB 8.6.C goldens + manifest CSV rows
# (Charter Section 3.3, apertura Section 5.6-5.10, parcial D34
# canonization). Mirror of test-eb_goldens_8_6_B.R for the
# multivariate-parcial sub-phase.
#
# Structural tests (no Stan re-run): load each .rds and verify that
# (i) the file exists when bootstrapped, (ii) the recorded
# family / regime / path / status match the roster, (iii) the
# manifest CSV contains a matching row with sub_phase = "8.6.C",
# (iv) the recorded numerics are finite for fitable cases and the
# guard_class is one of the canonical error classes for guard cases.
# Bootstrap lives in data-raw/bootstrap_eb_goldens_C.R and is gated
# by GDPAR_BOOTSTRAP_EB_GOLDENS=1 + GDPAR_BOOTSTRAP_EB_GOLDENS_C=1.

.eb_golden_C_root <- function() {
  c2 <- file.path("tests", "testthat", "data")
  if (dir.exists(c2)) return(c2)
  file.path("data")
}

.load_eb_golden_C_rds <- function(family, path, regime) {
  fname <- .gdpar_golden_eb_rds_name_8_6_C(family, path, regime)
  rds <- file.path(.eb_golden_C_root(), fname)
  if (!file.exists(rds)) return(NULL)
  list(path = rds, data = readRDS(rds))
}

# ---- Section 1: roster coverage and naming -------------------------------

test_that(".gdpar_golden_roster_8_6_C returns exactly 17 entries", {
  expect_length(.gdpar_golden_roster_8_6_C(), 17L)
})

test_that(".gdpar_golden_roster_8_6_C has 13 fitable + 4 guard cases", {
  roster <- .gdpar_golden_roster_8_6_C()
  statuses <- vapply(roster, function(x) x$status, character(1L))
  expect_equal(sum(statuses == "fitable"), 13L)
  expect_equal(sum(statuses == "guard"), 4L)
})

test_that(".gdpar_golden_roster_8_6_C splits 10 Path A + 7 Path B", {
  roster <- .gdpar_golden_roster_8_6_C()
  paths <- vapply(roster, function(x) x$path, character(1L))
  expect_equal(sum(paths == "A"), 10L)
  expect_equal(sum(paths == "B"), 7L)
})

test_that(".gdpar_golden_eb_rds_name_8_6_C follows canonical convention", {
  expect_equal(
    .gdpar_golden_eb_rds_name_8_6_C("gaussian", "A", "polyP2"),
    "golden_eb_gaussian_pathA_polyP2.rds"
  )
  expect_equal(
    .gdpar_golden_eb_rds_name_8_6_C("beta", "B", "K2"),
    "golden_eb_beta_pathB_K2.rds"
  )
  expect_equal(
    .gdpar_golden_eb_rds_name_8_6_C("beta", "A", "polyP2_guard"),
    "golden_eb_beta_pathA_polyP2_guard.rds"
  )
})

test_that("EB 8.6.C roster guard_class is canonical for every guard entry", {
  roster <- .gdpar_golden_roster_8_6_C()
  allowed <- c("gdpar_unsupported_feature_error",
               "gdpar_input_error")
  for (cfg in roster) {
    if (cfg$status == "guard") {
      expect_true(cfg$guard_class %in% allowed,
                  info = sprintf("guard_class %s not canonical for %s/%s/%s",
                                 cfg$guard_class, cfg$path,
                                 cfg$family, cfg$regime))
    } else {
      expect_true(is.na(cfg$guard_class),
                  info = sprintf("fitable %s/%s/%s should have NA guard_class",
                                 cfg$path, cfg$family, cfg$regime))
    }
  }
})

test_that("EB 8.6.C roster Path A K=1 entries have p > 1 and Path B K>1 have p=1", {
  roster <- .gdpar_golden_roster_8_6_C()
  for (cfg in roster) {
    if (cfg$path == "A") {
      expect_equal(cfg$K, 1L,
                   info = sprintf("Path A %s/%s K=%d (expected 1)",
                                  cfg$family, cfg$regime, cfg$K))
      expect_true(cfg$p > 1L,
                  info = sprintf("Path A %s/%s p=%d (expected > 1)",
                                 cfg$family, cfg$regime, cfg$p))
    } else if (cfg$path == "B") {
      expect_true(cfg$K > 1L,
                  info = sprintf("Path B %s/%s K=%d (expected > 1)",
                                 cfg$family, cfg$regime, cfg$K))
      expect_equal(cfg$p, 1L,
                   info = sprintf("Path B %s/%s p=%d (expected 1)",
                                  cfg$family, cfg$regime, cfg$p))
    }
  }
})

# ---- Section 2: golden RDS file presence (skip if not bootstrapped) -----

.eb_C_first_golden_present <- function() {
  any(vapply(.gdpar_golden_roster_8_6_C(), function(cfg) {
    file.exists(file.path(.eb_golden_C_root(),
                          .gdpar_golden_eb_rds_name_8_6_C(cfg$family,
                                                          cfg$path,
                                                          cfg$regime)))
  }, logical(1L)))
}

test_that("all 17 EB 8.6.C golden RDS files exist on disk (skip if not bootstrapped)", {
  testthat::skip_if_not(.eb_C_first_golden_present(),
                        "no EB 8.6.C goldens present (bootstrap not yet run)")
  for (cfg in .gdpar_golden_roster_8_6_C()) {
    rds_path <- file.path(.eb_golden_C_root(),
                          .gdpar_golden_eb_rds_name_8_6_C(cfg$family,
                                                          cfg$path,
                                                          cfg$regime))
    expect_true(file.exists(rds_path),
                info = sprintf("missing golden: %s", rds_path))
  }
})

# ---- Section 3: per-config RDS structural integrity (when present) ------

test_that("fitable EB 8.6.C goldens carry finite theta_ref_hat / theta_ref_se", {
  testthat::skip_if_not(.eb_C_first_golden_present(),
                        "no EB 8.6.C goldens present")
  for (cfg in .gdpar_golden_roster_8_6_C()) {
    if (cfg$status != "fitable") next
    loaded <- .load_eb_golden_C_rds(cfg$family, cfg$path, cfg$regime)
    if (is.null(loaded)) {
      testthat::skip(sprintf("not bootstrapped: %s/%s/%s",
                             cfg$path, cfg$family, cfg$regime))
    }
    d <- loaded$data
    expect_equal(d$family, cfg$family)
    expect_equal(d$path, cfg$path)
    expect_equal(d$regime, cfg$regime)
    expect_equal(d$status, "fitable")
    expect_equal(d$sub_phase, "8.6.C")
    expect_true(all(is.finite(d$theta_ref_hat)),
                info = sprintf("non-finite theta_ref_hat in %s/%s/%s",
                               cfg$path, cfg$family, cfg$regime))
    expect_true(all(is.finite(d$theta_ref_se)),
                info = sprintf("non-finite theta_ref_se in %s/%s/%s",
                               cfg$path, cfg$family, cfg$regime))
  }
})

test_that("guard EB 8.6.C goldens carry canonical guard_class metadata", {
  testthat::skip_if_not(.eb_C_first_golden_present(),
                        "no EB 8.6.C goldens present")
  for (cfg in .gdpar_golden_roster_8_6_C()) {
    if (cfg$status != "guard") next
    loaded <- .load_eb_golden_C_rds(cfg$family, cfg$path, cfg$regime)
    if (is.null(loaded)) next
    d <- loaded$data
    expect_equal(d$status, "guard")
    expect_equal(d$family, cfg$family)
    expect_equal(d$path, cfg$path)
    expect_equal(d$guard_class, cfg$guard_class)
  }
})

# ---- Section 4: manifest CSV cross-validation ----------------------------

test_that("manifest CSV contains sub_phase = '8.6.C' rows when bootstrapped", {
  testthat::skip_if_not(file.exists(.gdpar_golden_manifest_path()),
                        "golden_manifest.csv not present")
  testthat::skip_if_not(.eb_C_first_golden_present(),
                        "no EB 8.6.C goldens present")
  m <- .gdpar_golden_manifest_read()
  rows_C <- m[m$sub_phase == "8.6.C", , drop = FALSE]
  expect_true(nrow(rows_C) >= 1L)
})

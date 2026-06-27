# Sub-phase 8.6.D (Session 13c closure, 2026-05-26): structural
# validation of the Path C EB goldens + manifest CSV rows
# (Charter Section 3.4, apertura + parcial Session 13b, D36 = alpha +
# D37 = (i) + D38'' = (h) + D39 + D40' + D41 + D43 = (a) canonized).
# Mirror of test-eb_goldens_8_6_C.R for the K > 1 + p > 1 regime.
#
# Structural tests (no Stan re-run): load each .rds and verify that
# (i) the file exists when bootstrapped, (ii) the recorded family /
# regime / path / status / K / p match the roster, (iii) the manifest
# CSV contains a matching row with sub_phase = "8.6.D", (iv) the
# recorded numerics are finite for fitable cases and the guard_class
# is canonical for guarded ones. Bootstrap lives in
# data-raw/bootstrap_eb_goldens_D.R and is gated by
# GDPAR_BOOTSTRAP_EB_GOLDENS=1 + GDPAR_BOOTSTRAP_EB_GOLDENS_D=1.

.eb_golden_D_root <- function() {
  c2 <- file.path("tests", "testthat", "data")
  if (dir.exists(c2)) return(c2)
  file.path("data")
}

.load_eb_golden_D_rds <- function(family, regime) {
  fname <- .gdpar_golden_eb_rds_name_8_6_D(family, regime)
  rds <- file.path(.eb_golden_D_root(), fname)
  if (!file.exists(rds)) return(NULL)
  list(path = rds, data = readRDS(rds))
}

# ---- Section 1: roster coverage and naming -------------------------------

test_that(".gdpar_golden_roster_8_6_D returns exactly 8 entries", {
  expect_length(.gdpar_golden_roster_8_6_D(), 8L)
})

test_that(".gdpar_golden_roster_8_6_D has 4 fitable + 4 guard cases", {
  roster <- .gdpar_golden_roster_8_6_D()
  statuses <- vapply(roster, function(x) x$status, character(1L))
  expect_equal(sum(statuses == "fitable"), 4L)
  expect_equal(sum(statuses == "guard"), 4L)
})

test_that(".gdpar_golden_roster_8_6_D is all Path C", {
  roster <- .gdpar_golden_roster_8_6_D()
  paths <- vapply(roster, function(x) x$path, character(1L))
  expect_true(all(paths == "C"))
})

test_that(".gdpar_golden_roster_8_6_D enforces K > 1 AND p > 1 in every entry", {
  roster <- .gdpar_golden_roster_8_6_D()
  for (cfg in roster) {
    expect_true(cfg$K > 1L,
                info = sprintf("Path C %s/%s K=%d (expected > 1)",
                               cfg$family, cfg$regime, cfg$K))
    expect_true(cfg$p > 1L,
                info = sprintf("Path C %s/%s p=%d (expected > 1)",
                               cfg$family, cfg$regime, cfg$p))
  }
})

test_that(".gdpar_golden_roster_8_6_D fitable entries cover Gaussian + NB across p in {2, 3}", {
  roster <- .gdpar_golden_roster_8_6_D()
  fitable <- Filter(function(cfg) cfg$status == "fitable", roster)
  fams <- vapply(fitable, function(x) x$family, character(1L))
  ps <- vapply(fitable, function(x) x$p, integer(1L))
  expect_setequal(unique(fams), c("gaussian", "neg_binomial_2"))
  expect_setequal(unique(ps), c(2L, 3L))
})

test_that(".gdpar_golden_eb_rds_name_8_6_D follows canonical convention", {
  expect_equal(
    .gdpar_golden_eb_rds_name_8_6_D("gaussian", "polyP2"),
    "golden_eb_gaussian_pathC_polyP2.rds"
  )
  expect_equal(
    .gdpar_golden_eb_rds_name_8_6_D("neg_binomial_2", "polyP3"),
    "golden_eb_neg_binomial_2_pathC_polyP3.rds"
  )
  expect_equal(
    .gdpar_golden_eb_rds_name_8_6_D("beta", "polyP2_guard"),
    "golden_eb_beta_pathC_polyP2_guard.rds"
  )
})

test_that("EB 8.6.D roster guard_class is canonical for every guard entry", {
  roster <- .gdpar_golden_roster_8_6_D()
  allowed <- c("gdpar_unsupported_feature_error",
               "gdpar_input_error")
  for (cfg in roster) {
    if (cfg$status == "guard") {
      expect_true(cfg$guard_class %in% allowed,
                  info = sprintf("guard_class %s not canonical for C/%s/%s",
                                 cfg$guard_class, cfg$family, cfg$regime))
    } else {
      expect_true(is.na(cfg$guard_class),
                  info = sprintf("fitable C/%s/%s should have NA guard_class",
                                 cfg$family, cfg$regime))
    }
  }
})

# ---- Section 2: structural canonization invariants -----------------------

test_that(".gdpar_eb_correction_tensor shape is [K, p, p] when applied", {
  slots <- list(
    list(theta_ref_cov_k = matrix(c(0.5, 0.1, 0.1, 0.4), 2L, 2L)),
    list(theta_ref_cov_k = matrix(c(0.3, 0.0, 0.0, 0.2), 2L, 2L))
  )
  res <- .gdpar_eb_correction_tensor(
    eb_correction = TRUE,
    laplace_result_per_slot = slots,
    K = 2L, p = 2L, verbose = FALSE
  )
  expect_true(res$applied)
  expect_equal(dim(res$constant), c(2L, 2L, 2L))
  expect_equal(unname(res$slot_dispositions), c("ok", "ok"))
})

test_that(".gdpar_eb_correction_tensor with non-PSD slot yields applied = FALSE", {
  slots <- list(
    list(theta_ref_cov_k = matrix(c(0.5, 0.1, 0.1, 0.4), 2L, 2L)),
    list(theta_ref_cov_k = matrix(c(-1.0, 0.0, 0.0, -1.0), 2L, 2L))
  )
  res <- suppressWarnings(.gdpar_eb_correction_tensor(
    eb_correction = TRUE,
    laplace_result_per_slot = slots,
    K = 2L, p = 2L, verbose = FALSE
  ))
  expect_false(res$applied)
  expect_true("non_psd" %in% res$slot_dispositions)
})

test_that(".gdpar_eb_correction_tensor with eb_correction = FALSE yields disabled dispositions", {
  slots <- list(list(theta_ref_cov_k = matrix(0.5, 1L, 1L)),
                list(theta_ref_cov_k = matrix(0.3, 1L, 1L)))
  res <- .gdpar_eb_correction_tensor(
    eb_correction = FALSE,
    laplace_result_per_slot = slots,
    K = 2L, p = 1L, verbose = FALSE
  )
  expect_false(res$applied)
  expect_true(all(res$slot_dispositions == "disabled"))
})

test_that(".gdpar_eb_check_stan_id_for_path Path C accepts {1, 3} and rejects {5, 8}", {
  fam_g <- gdpar_family("gaussian")
  fam_n <- gdpar_family("neg_binomial_2")
  fam_b <- gdpar_family("beta")
  fam_s <- gdpar_family("student_t")
  expect_invisible(.gdpar_eb_check_stan_id_for_path(fam_g, 2L, 2L))
  expect_invisible(.gdpar_eb_check_stan_id_for_path(fam_n, 2L, 2L))
  expect_error(.gdpar_eb_check_stan_id_for_path(fam_b, 2L, 2L),
               class = "gdpar_unsupported_feature_error")
  expect_error(.gdpar_eb_check_stan_id_for_path(fam_s, 2L, 2L),
               class = "gdpar_unsupported_feature_error")
})

test_that(".gdpar_eb_check_stan_id_for_path Path B still accepts {5, 6, 7, 8, 9, 10, 11, 12, 13} at p = 1", {
  for (fam_name in c("beta", "gamma", "student_t", "tweedie",
                      "zip", "hurdle_poisson", "hurdle_neg_binomial_2")) {
    fam <- gdpar_family(fam_name)
    expect_invisible(.gdpar_eb_check_stan_id_for_path(fam, 2L, 1L),
                     label = paste0("Path B accept ", fam_name))
  }
})

# ---- Section 3: golden RDS file presence (skip if not bootstrapped) -----

.eb_D_first_golden_present <- function() {
  any(vapply(.gdpar_golden_roster_8_6_D(), function(cfg) {
    file.exists(file.path(.eb_golden_D_root(),
                          .gdpar_golden_eb_rds_name_8_6_D(cfg$family,
                                                          cfg$regime)))
  }, logical(1L)))
}

test_that("all EB 8.6.D golden RDS files exist on disk (skip if not bootstrapped)", {
  testthat::skip_if_not(.eb_D_first_golden_present(),
                        "no EB 8.6.D goldens present (bootstrap not yet run)")
  for (cfg in .gdpar_golden_roster_8_6_D()) {
    rds_path <- file.path(.eb_golden_D_root(),
                          .gdpar_golden_eb_rds_name_8_6_D(cfg$family,
                                                          cfg$regime))
    expect_true(file.exists(rds_path),
                info = sprintf("missing golden: %s", rds_path))
  }
})

# ---- Section 4: per-config RDS structural integrity (when present) ------

test_that("fitable EB 8.6.D goldens carry finite theta_ref_kp_hat and correction tensor", {
  testthat::skip_if_not(.eb_D_first_golden_present(),
                        "no EB 8.6.D goldens present")
  for (cfg in .gdpar_golden_roster_8_6_D()) {
    if (cfg$status != "fitable") next
    loaded <- .load_eb_golden_D_rds(cfg$family, cfg$regime)
    if (is.null(loaded)) {
      testthat::skip(sprintf("not bootstrapped: C/%s/%s",
                             cfg$family, cfg$regime))
    }
    d <- loaded$data
    expect_equal(d$family, cfg$family)
    expect_equal(d$path, "C")
    expect_equal(d$regime, cfg$regime)
    expect_equal(d$status, "fitable")
    expect_equal(d$sub_phase, "8.6.D")
    expect_equal(d$K, cfg$K)
    expect_equal(d$p, cfg$p)
    expect_true(all(is.finite(d$theta_ref_kp_hat)),
                info = sprintf("non-finite theta_ref_kp_hat in C/%s/%s",
                               cfg$family, cfg$regime))
    expect_equal(dim(d$correction_tensor_constant),
                 c(cfg$K, cfg$p, cfg$p))
  }
})

test_that("guard EB 8.6.D goldens carry canonical guard_class metadata", {
  testthat::skip_if_not(.eb_D_first_golden_present(),
                        "no EB 8.6.D goldens present")
  for (cfg in .gdpar_golden_roster_8_6_D()) {
    if (cfg$status != "guard") next
    loaded <- .load_eb_golden_D_rds(cfg$family, cfg$regime)
    if (is.null(loaded)) next
    d <- loaded$data
    expect_equal(d$status, "guard")
    expect_equal(d$family, cfg$family)
    expect_equal(d$path, "C")
    expect_equal(d$guard_class, cfg$guard_class)
  }
})

# ---- Section 5: manifest CSV cross-validation ----------------------------

test_that("manifest CSV contains sub_phase = '8.6.D' rows when bootstrapped", {
  testthat::skip_if_not(file.exists(.gdpar_golden_manifest_path()),
                        "golden_manifest.csv not present")
  testthat::skip_if_not(.eb_D_first_golden_present(),
                        "no EB 8.6.D goldens present")
  m <- .gdpar_golden_manifest_read()
  rows_D <- m[m$sub_phase == "8.6.D", , drop = FALSE]
  expect_true(nrow(rows_D) >= 1L)
})

# Sub-phase 8.3.10 (decision E6.C): bit-exact compare-path tests for
# the 14 goldens bootstrapped in 8.3.9. Tiered by computational cost
# so CI normal vs CI nocturno can opt into a partial vs full run:
#
#   Tier1 (~9 min, 12 configs, individual time <= 60s each except cfg 12
#   bspline_K2 at ~140s):
#     env var GDPAR_GOLDEN_8_3_9_COMPARE_TIER1=1 enables.
#   Tier2 (~40 min, 2 configs, Tweedie K=3 ~38min and bspline_p2 ~125s):
#     env var GDPAR_GOLDEN_8_3_9_COMPARE_TIER2=1 enables.
#
# Each test re-runs the original cfg_xx closure (sourced from
# inst/extdata/scripts/bootstrap_8_3_9_goldens.R, which now skips its
# top-level driver when source()-ed without GDPAR_BOOTSTRAP_8_3_9), then
# compares the resulting posterior draws_matrix bit-exactly against the
# .rds frozen during 8.3.9.
#
# Reproducibility: closures fix their own set.seed() and the bootstrap
# script pins CHAINS/ITER_WARMUP/ITER_SAMPLING. Bit-exactness within a
# fixed cmdstan version is the contract; cross-version drift in cmdstan
# 2.38.x or beyond is the indicator that the goldens need a deliberate
# re-bootstrap (gotcha 110 of handoff 38).

.skip_if_no_compare_tier <- function(tier) {
  envvar <- paste0("GDPAR_GOLDEN_8_3_9_COMPARE_TIER", tier)
  if (!nzchar(Sys.getenv(envvar))) {
    skip(sprintf("%s not set; skipping Tier%d compare-path", envvar, tier))
  }
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
}

.bootstrap_script_path <- function() {
  # Look for the source-tree path first (pkgload::load_all() workflow)
  # and fall back to the installed package layout (system.file()).
  # normalizePath() resolves relative segments so that the resulting
  # path is robust to subsequent setwd() calls.
  source_tree <- suppressWarnings(normalizePath(
    testthat::test_path("..", "..", "inst", "extdata", "scripts",
                        "bootstrap_8_3_9_goldens.R"),
    mustWork = FALSE
  ))
  candidates <- c(
    source_tree,
    system.file("extdata", "scripts", "bootstrap_8_3_9_goldens.R",
                package = "gdpar")
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(p)) return(p)
  }
  NA_character_
}

.load_bootstrap_closures <- function() {
  path <- .bootstrap_script_path()
  if (is.na(path)) {
    skip("bootstrap_8_3_9_goldens.R not found; skipping compare-path")
  }
  env <- new.env(parent = globalenv())
  # source() with skip_if_no banner: the script's top-level driver is
  # gated by .is_bootstrap_main_run() which returns FALSE when source()
  # is invoked from sys.nframe() > 0L (this test file). Run with wd set
  # to the package root so that any relative path in the script (e.g.,
  # cmdstan compilation cache, system.file lookups) resolves the way
  # the original bootstrap command-line invocation expects.
  pkg_root <- normalizePath(testthat::test_path("..", ".."),
                            mustWork = TRUE)
  old_wd <- setwd(pkg_root)
  on.exit(setwd(old_wd), add = TRUE)
  sys.source(path, envir = env)
  env
}

.compare_one_config <- function(label, env) {
  config_record <- NULL
  for (cfg in env$configs) {
    if (identical(cfg[[1L]], label)) {
      config_record <- cfg
      break
    }
  }
  expect_false(is.null(config_record),
               info = sprintf("config %s not found in bootstrap configs",
                              label))
  family     <- config_record[[3L]]
  K          <- config_record[[4L]]
  basis_type <- config_record[[5L]]
  p          <- config_record[[6L]]
  rds_name <- .gdpar_golden_rds_name(family, K, basis_type, p)
  rds_path <- testthat::test_path("data", rds_name)
  if (!file.exists(rds_path)) {
    skip(sprintf("Golden RDS missing for compare-path: %s", rds_path))
  }
  golden <- readRDS(rds_path)
  pkg_root <- normalizePath(testthat::test_path("..", ".."),
                            mustWork = TRUE)
  old_wd <- setwd(pkg_root)
  on.exit(setwd(old_wd), add = TRUE)
  fit <- config_record[[2L]]()
  expect_false(is.null(fit),
               info = sprintf("closure %s returned NULL", label))
  fresh_draws <- suppressWarnings(
    posterior::as_draws_matrix(fit$fit$draws())
  )
  expect_identical(
    unname(unclass(fresh_draws)),
    unname(unclass(golden$draws_matrix)),
    info = sprintf("bit-exact compare-path failed for %s", label)
  )
  invisible(TRUE)
}

# ---- Tier 1: 12 fast configs ---------------------------------------------

.TIER1_LABELS <- c(
  "01_lognormal_K2",
  "02_student_t_K3",
  "04_zip_K2",
  "05_zinb_K3",
  "06_hurdle_p_K2",
  "07_hurdle_nb_K3",
  "08_het_gauss_beta_K2",
  "09_het_gauss_gamma_K2",
  "10_het_nb_beta_K2",
  "11_gaussian_K1_bspline_p1",
  "12_het_gauss_beta_bspline_K2",
  "14_gaussian_K1_polynomial_p2"
)

test_that("Tier1 bit-exact compare-path: 12 fast configs match their frozen RDS draws (~9 min)", {
  .skip_if_no_compare_tier(1L)
  env <- .load_bootstrap_closures()
  for (label in .TIER1_LABELS) {
    .compare_one_config(label, env)
  }
})

# ---- Tier 2: 2 slow configs ---------------------------------------------

.TIER2_LABELS <- c(
  "03_tweedie_K3",           # ~38 min Stan sampling
  "13_gaussian_K1_bspline_p2" # ~125s Stan sampling
)

test_that("Tier2 bit-exact compare-path: 2 slow configs match their frozen RDS draws (~40 min)", {
  .skip_if_no_compare_tier(2L)
  env <- .load_bootstrap_closures()
  for (label in .TIER2_LABELS) {
    .compare_one_config(label, env)
  }
})

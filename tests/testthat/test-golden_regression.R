## Env-gated golden regression tests for Phase I.3.
##
## Each test refits the canonical smoke_p2_* scenario with the same
## seed and data as inst/benchmarks/scripts/build_golden_multi.R,
## extracts the observed snapshot via gdpar_snapshot_fit(),
## and compares it against the persisted golden .rds using the
## four-layer comparator gdpar_golden_compare(). All four layers
## are evaluated and the failure rows are surfaced through testthat's
## info string on failure.
##
## Wall time: about 5-7 minutes per test (one Stan compile cached
## across tests; three short fits of 1000/1000). Use the env gate
## GDPAR_GOLDEN_CHECK=1 to run them.

skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  has_cs <- tryCatch(
    {
      cmdstanr::cmdstan_version()
      TRUE
    },
    error = function(e) FALSE
  )
  testthat::skip_if(!has_cs, "cmdstan not available")
}

build_smoke_p2_data_test <- function(seed = 42L, n = 300L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  beta_a_true <- list(c(2.0, -1.0), c(-2.0, 1.0))
  theta_ref_true <- c(0.5, -0.5)
  sigma_true <- 0.3
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1_c +
      beta_a_true[[1]][2] * x2_c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1_c +
      beta_a_true[[2]][2] * x2_c
  )
  y_mat <- eta_true + matrix(rnorm(n * 2L, 0, sigma_true), n, 2L)
  df <- data.frame(x1 = x1, x2 = x2)
  df$y <- y_mat
  df
}

run_golden_regression_case <- function(scenario_id, parametrization) {
  golden_path <- system.file(
    "benchmarks", "golden", paste0(scenario_id, ".rds"),
    package = "gdpar"
  )
  if (!nzchar(golden_path) || !file.exists(golden_path)) {
    golden_path <- file.path(
      "inst", "benchmarks", "golden",
      paste0(scenario_id, ".rds")
    )
  }
  testthat::skip_if(
    !file.exists(golden_path),
    sprintf(paste(
      "Golden not present at %s.",
      "Run Rscript inst/benchmarks/scripts/build_golden_multi.R to",
      "generate the goldens first."
    ), golden_path)
  )
  golden <- readRDS(golden_path)
  df <- build_smoke_p2_data_test(seed = golden$seed, n = golden$n)
  spec <- amm_spec(p = 2L,
                    dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  fit <- suppressMessages(suppressWarnings(
    gdpar(
      formula       = y ~ x1 + x2,
      family        = fam, amm = spec, data = df,
      chains        = golden$n_chains,
      iter_warmup   = golden$n_warmup,
      iter_sampling = golden$n_sampling,
      refresh       = 0L, seed = golden$seed,
      parametrization = parametrization,
      skip_id_check = FALSE, verbose = FALSE
    )
  ))
  snapshot <- gdpar_snapshot_fit(fit)
  cmp <- gdpar_golden_compare(snapshot, golden, k_sigma = 3)
  fail_info <- if (!cmp$passed) {
    paste0(
      "\nFailures by layer: ",
      paste(names(cmp$by_layer), cmp$by_layer, sep = "=",
             collapse = ", "),
      "\nFailure rows:\n",
      paste(utils::capture.output(print(cmp$failures,
                                          row.names = FALSE)),
            collapse = "\n")
    )
  } else {
    ""
  }
  testthat::expect_true(cmp$passed, info = fail_info)
}

test_that("golden_regression: smoke_p2_auto", {
  testthat::skip_if(Sys.getenv("GDPAR_GOLDEN_CHECK") != "1",
    "Set GDPAR_GOLDEN_CHECK=1 to run the multi golden regression (compiles a Stan model; takes minutes).")
  skip_if_no_cmdstan()
  run_golden_regression_case("smoke_p2_auto", "auto")
})

test_that("golden_regression: smoke_p2_cp", {
  testthat::skip_if(Sys.getenv("GDPAR_GOLDEN_CHECK") != "1",
    "Set GDPAR_GOLDEN_CHECK=1 to run the multi golden regression (compiles a Stan model; takes minutes).")
  skip_if_no_cmdstan()
  run_golden_regression_case("smoke_p2_cp", "cp")
})

test_that("golden_regression: smoke_p2_ncp", {
  testthat::skip_if(Sys.getenv("GDPAR_GOLDEN_CHECK") != "1",
    "Set GDPAR_GOLDEN_CHECK=1 to run the multi golden regression (compiles a Stan model; takes minutes).")
  skip_if_no_cmdstan()
  run_golden_regression_case("smoke_p2_ncp", "ncp")
})

# ---- pure tests of gdpar_golden_compare layers (always run) ----------

test_that("gdpar_golden_compare passes when observed equals golden", {
  golden <- list(
    structural = list(
      fit_class = "gdpar_fit", coef_class = "gdpar_coef",
      p = 2L, summary_stats = c("mean", "q05"),
      components = c("theta_ref", "a", "b", "W"),
      theta_ref_cols = c("k", "mean", "q05"),
      theta_ref_nrow = 2L,
      a_class = "list", a_length = 2L,
      a_per_k_cols = list(c("term", "mean"), c("term", "mean")),
      b_class = "NULL", b_length = NA_integer_,
      W_class = "NULL", W_length = NA_integer_,
      parametrization_keys = c("cp_a", "cp_W")
    ),
    discrete = list(
      n_divergent = 0L, treedepth_max_n = 0L,
      treedepth_max_value = 7L,
      n_leapfrog_total_per_chain = c(1000L, 1100L),
      ebfmi_min = 0.8
    ),
    continuous = list(
      theta_ref = list(
        "theta_ref[1]" = list(
          mean = 0.5, sd = 0.01, ess_bulk = 500,
          ess_tail = 600, rhat = 1.001, mc_se = 0.0005
        )
      )
    ),
    sanity = list(
      rhat_max = 1.005, ess_bulk_min = 400,
      ess_tail_min = 500, divergent_pct = 0, ebfmi_min = 0.8
    )
  )
  cmp <- gdpar_golden_compare(golden, golden)
  expect_true(cmp$passed)
  expect_equal(nrow(cmp$failures), 0L)
})

test_that("gdpar_golden_compare flags structural drift", {
  golden <- list(
    structural = list(
      fit_class = "gdpar_fit", coef_class = "gdpar_coef",
      p = 2L, summary_stats = c("mean", "q05"),
      components = c("theta_ref", "a", "b", "W"),
      theta_ref_cols = c("k", "mean", "q05"),
      theta_ref_nrow = 2L,
      a_class = "list", a_length = 2L,
      a_per_k_cols = list(c("term", "mean"), c("term", "mean")),
      b_class = "NULL", b_length = NA_integer_,
      W_class = "NULL", W_length = NA_integer_,
      parametrization_keys = c("cp_a", "cp_W")
    ),
    discrete = list(
      n_divergent = 0L, treedepth_max_n = 0L,
      treedepth_max_value = 7L,
      n_leapfrog_total_per_chain = c(1000L, 1100L),
      ebfmi_min = 0.8
    ),
    continuous = list(),
    sanity = list(
      rhat_max = 1.005, ess_bulk_min = 400,
      ess_tail_min = 500, divergent_pct = 0, ebfmi_min = 0.8
    )
  )
  observed <- golden
  observed$structural$p <- 3L
  cmp <- gdpar_golden_compare(observed, golden)
  expect_false(cmp$passed)
  expect_true(any(cmp$failures$layer == "structural"))
  expect_true(any(cmp$failures$item == "p"))
})

test_that("gdpar_golden_compare flags continuous drift beyond k_sigma * MC_SE", {
  cont_golden <- list(theta_ref = list(
    "theta_ref[1]" = list(mean = 0.5, sd = 0.01, ess_bulk = 500,
                           ess_tail = 600, rhat = 1.001,
                           mc_se = 0.001)
  ))
  cont_obs <- list(theta_ref = list(
    "theta_ref[1]" = list(mean = 0.51, sd = 0.01, ess_bulk = 500,
                           ess_tail = 600, rhat = 1.001,
                           mc_se = 0.001)
  ))
  golden_min <- list(
    structural = NULL, discrete = NULL,
    continuous = cont_golden,
    sanity = list(rhat_max = 1.0, ess_bulk_min = 1000,
                   ess_tail_min = 1000, divergent_pct = 0,
                   ebfmi_min = 1.0)
  )
  obs_min <- list(
    structural = NULL, discrete = NULL,
    continuous = cont_obs,
    sanity = list(rhat_max = 1.0, ess_bulk_min = 1000,
                   ess_tail_min = 1000, divergent_pct = 0,
                   ebfmi_min = 1.0)
  )
  cmp <- gdpar_golden_compare(obs_min, golden_min, k_sigma = 3)
  expect_false(cmp$passed)
  expect_true(any(cmp$failures$layer == "continuous"))
})

test_that("gdpar_golden_compare flags sanity floor violations", {
  golden_min <- list(
    structural = NULL, discrete = NULL, continuous = NULL,
    sanity = list(rhat_max = 1.0, ess_bulk_min = 1000,
                   ess_tail_min = 1000, divergent_pct = 0,
                   ebfmi_min = 1.0)
  )
  obs_min <- list(
    structural = NULL, discrete = NULL, continuous = NULL,
    sanity = list(rhat_max = 1.10, ess_bulk_min = 50,
                   ess_tail_min = 50, divergent_pct = 0.05,
                   ebfmi_min = 0.1)
  )
  cmp <- gdpar_golden_compare(obs_min, golden_min)
  expect_false(cmp$passed)
  expect_true(any(cmp$failures$layer == "sanity"))
  expect_true(all(c("rhat_max", "ess_bulk_min", "ess_tail_min",
                     "divergent_pct", "ebfmi_min") %in%
                   cmp$failures$item))
})

test_that("gdpar_golden_compare flags discrete bit-exact mismatch", {
  disc_golden <- list(
    n_divergent = 0L, treedepth_max_n = 0L,
    treedepth_max_value = 7L,
    n_leapfrog_total_per_chain = c(1000L, 1100L),
    ebfmi_min = 0.8
  )
  disc_obs <- disc_golden
  disc_obs$n_divergent <- 3L
  golden_min <- list(
    structural = NULL, discrete = disc_golden, continuous = NULL,
    sanity = list(rhat_max = 1.0, ess_bulk_min = 1000,
                   ess_tail_min = 1000, divergent_pct = 0,
                   ebfmi_min = 1.0)
  )
  obs_min <- list(
    structural = NULL, discrete = disc_obs, continuous = NULL,
    sanity = list(rhat_max = 1.0, ess_bulk_min = 1000,
                   ess_tail_min = 1000, divergent_pct = 0,
                   ebfmi_min = 1.0)
  )
  cmp <- gdpar_golden_compare(obs_min, golden_min)
  expect_false(cmp$passed)
  expect_true(any(cmp$failures$layer == "discrete"))
  expect_true("n_divergent" %in% cmp$failures$item)
})

test_that("format_compact handles all input flavors", {
  expect_equal(gdpar:::format_compact(NULL), "<NULL>")
  expect_equal(gdpar:::format_compact(integer(0)), "<empty>")
  expect_equal(gdpar:::format_compact(c(1.5, 2.5)), "1.5,2.5")
  expect_equal(gdpar:::format_compact(c("a", "b")), "a,b")
  expect_equal(gdpar:::format_compact(c(TRUE, FALSE)), "TRUE,FALSE")
  expect_true(grepl("<num\\[", gdpar:::format_compact(rnorm(10))))
})

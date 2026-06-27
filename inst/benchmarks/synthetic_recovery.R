# ============================================================
# gdpar / inst/benchmarks/synthetic_recovery.R
# Parameter-recovery tests on synthetic data with known truth.
# Verifies that gdpar recovers the data-generating parameters
# (theta_ref, basis coefficients) within Monte Carlo uncertainty.
#
# Designed to take ~10-20 minutes total on a modern multi-core CPU.
# Set REPLICATES higher for tighter empirical bounds at higher cost.
# ============================================================

if (!exists("bench_section", mode = "function")) {
  helpers_path <- local({
    candidates <- character(0)
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      ctx <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                      error = function(e) "")
      if (nzchar(ctx) && file.exists(ctx)) {
        candidates <- c(candidates,
                        file.path(dirname(ctx), "helpers.R"))
      }
    }
    sourced <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(sourced) && nzchar(sourced)) {
      candidates <- c(candidates,
                      file.path(dirname(normalizePath(sourced,
                                                     mustWork = FALSE)),
                                "helpers.R"))
    }
    candidates <- c(candidates,
                    file.path(getwd(), "helpers.R"),
                    file.path(getwd(), "inst", "benchmarks",
                              "helpers.R"))
    inst_path <- system.file("benchmarks", "helpers.R",
                             package = "gdpar")
    if (nzchar(inst_path)) candidates <- c(candidates, inst_path)
    found <- ""
    for (p in candidates) {
      if (nzchar(p) && file.exists(p)) { found <- p; break }
    }
    if (!nzchar(found)) {
      stop("Could not locate helpers.R. Set the working directory ",
           "to the gdpar package root or to its inst/benchmarks ",
           "subdirectory, or install gdpar so that system.file() ",
           "resolves it.")
    }
    found
  })
  source(helpers_path)
}

bench_section("gdpar synthetic parameter-recovery tests")

if (!bench_require("gdpar")) stop("gdpar must be installed first.")
library(gdpar)
have_cmdstan <- bench_have_cmdstan()
if (!have_cmdstan) {
  bench_skip("All recovery tests require cmdstan; aborting.")
  invisible(NULL)
} else {

REPLICATES <- 5L
SAMPLE_SIZES <- c(200L, 500L, 1000L)
RESULTS <- list()

# ------------------------------------------------------------
# Recovery 1: Gaussian Level 1 -- recovery of theta_ref and a
# ------------------------------------------------------------

bench_subsection("Recovery 1: Gaussian Level 1, varying n")

for (n in SAMPLE_SIZES) {
  bench_status(sprintf("n = %d", n))
  rec <- vector("list", REPLICATES)
  for (r in seq_len(REPLICATES)) {
    df <- make_gauss_clean(n = n, seed = 1000L + r)
    fit <- gdpar(
      formula       = y ~ x1 + x2,
      family        = gdpar_family("gaussian"),
      amm           = amm_spec(a = ~ x1 + x2),
      data          = df,
      iter_warmup   = 500,
      iter_sampling = 500,
      chains        = 2,
      refresh       = 0,
      verbose       = FALSE
    )
    cf <- coef(fit)
    rec[[r]] <- list(
      theta_ref = cf$theta_ref,
      a_x1 = cf$a["x1"],
      a_x2 = cf$a["x2"],
      rhat_max = fit$diagnostics$rhat_max,
      ess_min = fit$diagnostics$ess_bulk_min,
      diverg = fit$diagnostics$divergent_count
    )
  }
  est <- do.call(rbind, lapply(rec, as.data.frame))
  truth <- c(theta_ref = 1, a_x1 = 0.4, a_x2 = -0.4)
  bias <- c(
    theta_ref = mean(est$theta_ref) - unname(truth["theta_ref"]),
    a_x1 = mean(est$a_x1) - unname(truth["a_x1"]),
    a_x2 = mean(est$a_x2) - unname(truth["a_x2"])
  )
  rmse_recovery <- c(
    theta_ref = sqrt(mean((est$theta_ref - unname(truth["theta_ref"]))^2)),
    a_x1 = sqrt(mean((est$a_x1 - unname(truth["a_x1"]))^2)),
    a_x2 = sqrt(mean((est$a_x2 - unname(truth["a_x2"]))^2))
  )
  RESULTS[[paste0("gauss_l1_n", n)]] <- list(
    n = n,
    bias = bias,
    rmse = rmse_recovery,
    mean_rhat_max = mean(est$rhat_max),
    mean_ess_min = mean(est$ess_min),
    total_divergent = sum(est$diverg)
  )
  bench_status(sprintf(
    "  bias theta_ref = %+.4f, RMSE = %.4f; bias a_x1 = %+.4f, RMSE = %.4f",
    bias["theta_ref"], rmse_recovery["theta_ref"],
    bias["a_x1"], rmse_recovery["a_x1"]
  ))
}

# ------------------------------------------------------------
# Recovery 2: Gaussian Level 2 -- recovery of all components
# ------------------------------------------------------------

bench_subsection("Recovery 2: Gaussian Level 2 with active b, n = 1000")

theta_ref_true <- 1.5
a_true <- c(0.4, -0.4)
b_true <- c(0.3, -0.3)
W_true <- c(0, 0)  # W block disabled to keep recovery clean for a + b

rec_l2 <- vector("list", REPLICATES)
for (r in seq_len(REPLICATES)) {
  df <- make_gauss_amm_level2(n = 1000, theta_ref = theta_ref_true,
                               a_coef = a_true, b_coef = b_true,
                               W_coef = W_true,
                               seed = 2000L + r)
  fit <- gdpar(
    formula       = y ~ x1 + x2 + x3 + x4,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1 + x2, b = ~ x3 + x4),
    data          = df,
    iter_warmup   = 2000,
    iter_sampling = 1000,
    chains        = 4,
    adapt_delta   = 0.99,
    max_treedepth = 15,
    refresh       = 0,
    verbose       = FALSE,
    anchor        = "empirical_y"
  )
  cf <- coef(fit)
  bench_status(sprintf(
    "  rep %d: rhat=%.3f div=%d theta_ref=%.3f a_x1=%.3f a_x2=%.3f b_x3=%.3f b_x4=%.3f",
    r, fit$diagnostics$rhat_max, fit$diagnostics$divergent_count,
    cf$theta_ref, cf$a["x1"], cf$a["x2"], cf$b["x3"], cf$b["x4"]
  ))
  rec_l2[[r]] <- list(
    theta_ref = cf$theta_ref,
    a_x1 = cf$a["x1"], a_x2 = cf$a["x2"],
    b_x3 = cf$b["x3"], b_x4 = cf$b["x4"],
    rhat_max = fit$diagnostics$rhat_max,
    diverg = fit$diagnostics$divergent_count
  )
}
est_l2 <- do.call(rbind, lapply(rec_l2, as.data.frame))
bench_status("  Mean estimates over replicates:")
bench_status(sprintf(
  "    theta_ref: est = %.3f (truth %.3f)",
  mean(est_l2$theta_ref), theta_ref_true
))
bench_status(sprintf(
  "    a_x1: est = %.3f (truth %.3f)",
  mean(est_l2$a_x1), a_true[1]
))
bench_status(sprintf(
  "    a_x2: est = %.3f (truth %.3f)",
  mean(est_l2$a_x2), a_true[2]
))
bench_status(sprintf(
  "    b_x3: est = %.3f (truth %.3f)",
  mean(est_l2$b_x3), b_true[1]
))
bench_status(sprintf(
  "    b_x4: est = %.3f (truth %.3f)",
  mean(est_l2$b_x4), b_true[2]
))
RESULTS$gauss_l2 <- list(
  est = est_l2,
  truth = c(theta_ref = theta_ref_true, a = a_true, b = b_true)
)

# ------------------------------------------------------------
# Recovery 3: Poisson Level 1
# ------------------------------------------------------------

bench_subsection("Recovery 3: Poisson Level 1, n = 500")

rec_pois <- vector("list", REPLICATES)
for (r in seq_len(REPLICATES)) {
  df <- make_poisson_clean(n = 500, seed = 3000L + r)
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family("poisson"),
    amm           = amm_spec(a = ~ x1 + x2),
    data          = df,
    iter_warmup   = 500,
    iter_sampling = 500,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  cf <- coef(fit)
  rec_pois[[r]] <- list(
    theta_ref = cf$theta_ref,
    a_x1 = cf$a["x1"],
    a_x2 = cf$a["x2"],
    rhat_max = fit$diagnostics$rhat_max
  )
}
est_pois <- do.call(rbind, lapply(rec_pois, as.data.frame))
truth_pois <- c(theta_ref = 1, a_x1 = 0.25, a_x2 = -0.25)
bench_status(sprintf(
  "  bias theta_ref = %+.4f; bias a_x1 = %+.4f; bias a_x2 = %+.4f",
  mean(est_pois$theta_ref) - truth_pois["theta_ref"],
  mean(est_pois$a_x1) - truth_pois["a_x1"],
  mean(est_pois$a_x2) - truth_pois["a_x2"]
))
RESULTS$pois_l1 <- list(est = est_pois, truth = truth_pois)

# ------------------------------------------------------------
# Recovery 4: Negative Binomial Level 1
# ------------------------------------------------------------

bench_subsection("Recovery 4: Negative Binomial Level 1, n = 500")

rec_nb <- vector("list", REPLICATES)
for (r in seq_len(REPLICATES)) {
  df <- make_negbin_clean(n = 500, phi = 5, seed = 4000L + r)
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family("neg_binomial_2"),
    amm           = amm_spec(a = ~ x1 + x2),
    data          = df,
    iter_warmup   = 500,
    iter_sampling = 500,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  cf <- coef(fit)
  rec_nb[[r]] <- list(
    theta_ref = cf$theta_ref,
    a_x1 = cf$a["x1"],
    a_x2 = cf$a["x2"],
    rhat_max = fit$diagnostics$rhat_max
  )
}
est_nb <- do.call(rbind, lapply(rec_nb, as.data.frame))
truth_nb <- c(theta_ref = 1, a_x1 = 0.25, a_x2 = -0.25)
bench_status(sprintf(
  "  bias theta_ref = %+.4f; bias a_x1 = %+.4f; bias a_x2 = %+.4f",
  mean(est_nb$theta_ref) - truth_nb["theta_ref"],
  mean(est_nb$a_x1) - truth_nb["a_x1"],
  mean(est_nb$a_x2) - truth_nb["a_x2"]
))
RESULTS$negbin_l1 <- list(est = est_nb, truth = truth_nb)

# ------------------------------------------------------------
# Recovery 5: Bernoulli Level 1
# ------------------------------------------------------------

bench_subsection("Recovery 5: Bernoulli Level 1, n = 800")

rec_be <- vector("list", REPLICATES)
for (r in seq_len(REPLICATES)) {
  df <- make_bernoulli_clean(n = 800, seed = 5000L + r)
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family("bernoulli"),
    amm           = amm_spec(a = ~ x1 + x2),
    data          = df,
    iter_warmup   = 500,
    iter_sampling = 500,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  cf <- coef(fit)
  rec_be[[r]] <- list(
    theta_ref = cf$theta_ref,
    a_x1 = cf$a["x1"],
    a_x2 = cf$a["x2"],
    rhat_max = fit$diagnostics$rhat_max
  )
}
est_be <- do.call(rbind, lapply(rec_be, as.data.frame))
truth_be <- c(theta_ref = 0, a_x1 = 0.55, a_x2 = -0.55)
bench_status(sprintf(
  "  bias theta_ref = %+.4f; bias a_x1 = %+.4f; bias a_x2 = %+.4f",
  mean(est_be$theta_ref) - truth_be["theta_ref"],
  mean(est_be$a_x1) - truth_be["a_x1"],
  mean(est_be$a_x2) - truth_be["a_x2"]
))
RESULTS$bern_l1 <- list(est = est_be, truth = truth_be)

# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

bench_section("Recovery test summary")

table_rows <- list()
for (n in SAMPLE_SIZES) {
  res <- RESULTS[[paste0("gauss_l1_n", n)]]
  table_rows[[length(table_rows) + 1L]] <- data.frame(
    family = "gaussian", level = 1, n = res$n,
    bias_theta_ref = res$bias["theta_ref"],
    rmse_theta_ref = res$rmse["theta_ref"],
    mean_rhat = res$mean_rhat_max,
    mean_ess = res$mean_ess_min,
    total_div = res$total_divergent,
    stringsAsFactors = FALSE
  )
}
table_rows[[length(table_rows) + 1L]] <- data.frame(
  family = "gaussian", level = 2, n = 1000,
  bias_theta_ref = mean(RESULTS$gauss_l2$est$theta_ref) - theta_ref_true,
  rmse_theta_ref = sqrt(mean((RESULTS$gauss_l2$est$theta_ref - theta_ref_true)^2)),
  mean_rhat = mean(RESULTS$gauss_l2$est$rhat_max),
  mean_ess = NA_real_,
  total_div = sum(RESULTS$gauss_l2$est$diverg),
  stringsAsFactors = FALSE
)
table_rows[[length(table_rows) + 1L]] <- data.frame(
  family = "poisson", level = 1, n = 500,
  bias_theta_ref = mean(RESULTS$pois_l1$est$theta_ref) - 1,
  rmse_theta_ref = sqrt(mean((RESULTS$pois_l1$est$theta_ref - 1)^2)),
  mean_rhat = mean(RESULTS$pois_l1$est$rhat_max),
  mean_ess = NA_real_,
  total_div = NA_integer_,
  stringsAsFactors = FALSE
)
table_rows[[length(table_rows) + 1L]] <- data.frame(
  family = "neg_binomial_2", level = 1, n = 500,
  bias_theta_ref = mean(RESULTS$negbin_l1$est$theta_ref) - 1,
  rmse_theta_ref = sqrt(mean((RESULTS$negbin_l1$est$theta_ref - 1)^2)),
  mean_rhat = mean(RESULTS$negbin_l1$est$rhat_max),
  mean_ess = NA_real_,
  total_div = NA_integer_,
  stringsAsFactors = FALSE
)
table_rows[[length(table_rows) + 1L]] <- data.frame(
  family = "bernoulli", level = 1, n = 800,
  bias_theta_ref = mean(RESULTS$bern_l1$est$theta_ref) - 0,
  rmse_theta_ref = sqrt(mean((RESULTS$bern_l1$est$theta_ref - 0)^2)),
  mean_rhat = mean(RESULTS$bern_l1$est$rhat_max),
  mean_ess = NA_real_,
  total_div = NA_integer_,
  stringsAsFactors = FALSE
)
summary_df <- do.call(rbind, table_rows)
bench_table(summary_df, digits = 4)

cat("\nNotes:\n")
cat("- bias and RMSE are over", REPLICATES, "replicates per row.\n")
cat("- mean_rhat is the mean of the maximum R-hat across replicates.\n")
cat("- A converged fit has mean_rhat <= 1.01.\n")
cat("- Increase REPLICATES at the top of this script for tighter Monte Carlo bounds.\n")

assign("synthetic_recovery_results", RESULTS, envir = .GlobalEnv)
invisible(NULL)
}

## S7 diagnostic — gdpar vs brms shrinkage analysis.
##
## Block 6 (handoff 17): S7 (p=3 heterogeneous gaussian) showed gdpar
## ELPD = -1466.5 vs ~-1377 for the other four competitors (delta ~88 pt,
## outside 2*se = 70). Sampling diagnostics were clean (0 divergences,
## rhat 1.00, ess > 2200). Hypothesis from the handoff: restrictive
## defaults on the hierarchical scale priors cause differential shrinkage
## when coordinate-wise effect magnitudes are very heterogeneous.
##
## This script fits gdpar and brms (set_rescor=FALSE) on S7 once, then
## extracts side-by-side estimates of:
##   - coordinate-wise beta coefficients (gdpar b_coef vs brms fixef)
##   - hierarchical scales (gdpar sigma_a, sigma_b, sigma_W, sigma_y)
##   - theta_ref posterior vs ground truth
##   - per-coordinate pointwise ELPD (where does the 88-pt gap live?)
##
## Optionally, if GDPAR_S7_REFIT_WIDER=1, fits gdpar a second time with
## a wider prior on sigma_b and reports the ELPD delta — testing whether
## the regression to -1378 is achievable with a prior change alone.
##
## Run:
##   cd <gdpar repo>
##   Rscript inst/benchmarks/diagnose_s7.R
##
## With wider-prior refit (~10 min extra):
##   GDPAR_S7_REFIT_WIDER=1 Rscript inst/benchmarks/diagnose_s7.R
##
## Output:
##   inst/benchmarks/results/s7_diagnostic.md  (structured summary)
##   stdout                                     (full diagnostic log)

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressMessages({
  devtools::load_all(".")
  source("inst/benchmarks/scripts/bench_multi_helpers.R")
  source("inst/benchmarks/scripts/bench_competitor_gdpar.R")
  source("inst/benchmarks/scripts/bench_competitor_brms.R")
  source("inst/benchmarks/synthetic_hard_multi.R")
})

results_dir <- "inst/benchmarks/results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
report_path <- file.path(results_dir, "s7_diagnostic.md")
report <- character()
push <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  report <<- c(report, msg)
}

log_sum_exp <- function(z) {
  m <- max(z)
  m + log(sum(exp(z - m)))
}

push("# S7 diagnostic report")
push("")
push("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
push("")

sc <- gen_S7_heterog_p3()
true_params <- sc$true_params

push("## Ground truth (S7)")
push("")
push("- theta_ref     : ", paste(true_params$theta_ref, collapse = ", "))
push("- beta_a matrix (rows = x1, x2; cols = coord 1..3):")
for (r in seq_len(nrow(true_params$beta_a))) {
  push("    [", c("x1", "x2")[r], "]  ",
       paste(formatC(true_params$beta_a[r, ], digits = 3,
                       format = "f", flag = "+"),
             collapse = "  "))
}
push("- sigma_y       : ", paste(true_params$sigma_y, collapse = ", "))
push("")

push("## Fitting gdpar with default priors ...")
push("")
t0 <- Sys.time()
fit_g <- bench_gdpar_fit(sc$train, sc$family_spec, sc$scenario_meta)
gdpar_wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
push("gdpar wall_seconds: ", round(gdpar_wall, 1))
push("")

push("## Fitting brms (set_rescor = FALSE) ...")
push("")
t1 <- Sys.time()
fit_b <- bench_brms_fit(sc$train, sc$family_spec, sc$scenario_meta)
brms_wall <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
push("brms wall_seconds: ", round(brms_wall, 1))
push("")

## -- gdpar coef summary --------------------------------------------------
push("## gdpar coef() summary")
push("")
cg <- coef(fit_g$fit)
df_cg <- as.data.frame(cg)
push("```")
capture.output(print(df_cg, row.names = FALSE)) |> (\(z) {
  for (ln in z) push(ln)
})()
push("```")
push("")

## -- brms fixef ----------------------------------------------------------
push("## brms fixef()")
push("")
bf <- brms::fixef(fit_b$fit)
push("```")
capture.output(print(round(bf, 4))) |> (\(z) {
  for (ln in z) push(ln)
})()
push("```")
push("")

## -- side-by-side beta comparison ----------------------------------------
push("## Side-by-side beta comparison (gdpar b_coef vs brms vs true)")
push("")
p <- as.integer(sc$scenario_meta$p)

## as.data.frame.gdpar_coef columns: component, k, identifier, x_name,
## mean, q05, q50, q95. The S7 AMM is `dimwise(a = ~ x1 + x2, b = NULL)`,
## so the linear-predictor coefficients of x1, x2 live in the `a`
## component, not `b`. We extract a-rows and split (Intercept) vs covars.
gdpar_a <- if ("component" %in% names(df_cg))
  df_cg[df_cg$component == "a", , drop = FALSE] else NULL
gdpar_b <- if ("component" %in% names(df_cg))
  df_cg[df_cg$component == "b", , drop = FALSE] else NULL
if (!is.null(gdpar_a) && nrow(gdpar_a) > 0L) {
  push("gdpar a rows (intercept + covariate effects on linear predictor):")
  push("```")
  capture.output(print(gdpar_a, row.names = FALSE)) |> (\(z) {
    for (ln in z) push(ln)
  })()
  push("```")
}
if (!is.null(gdpar_b) && nrow(gdpar_b) > 0L) {
  push("")
  push("gdpar b rows (W-scale b_coef):")
  push("```")
  capture.output(print(gdpar_b, row.names = FALSE)) |> (\(z) {
    for (ln in z) push(ln)
  })()
  push("```")
}
push("")

push("True betas vs gdpar posterior mean vs brms posterior mean (per coord, per covariate):")
push("")
push("```")
push(sprintf("%-6s %-4s %10s %10s %10s %10s %10s",
              "coord", "var", "true", "gdpar_mean", "brms_mean",
              "gdpar_dev_pct", "brms_dev_pct"))
for (k in seq_len(p)) {
  for (j in seq_along(c("x1", "x2"))) {
    var_name <- c("x1", "x2")[j]
    true_val <- true_params$beta_a[j, k]
    brms_name <- sprintf("y%d_%s", k, var_name)
    brms_mean <- if (brms_name %in% rownames(bf))
      bf[brms_name, "Estimate"] else NA_real_
    gdpar_mean <- NA_real_
    if (!is.null(gdpar_a) && nrow(gdpar_a) > 0L) {
      msk <- gdpar_a$k == k &
              (gdpar_a$x_name == var_name |
               gdpar_a$identifier == var_name)
      if (any(msk)) gdpar_mean <- gdpar_a$mean[msk][1L]
    }
    g_dev <- if (!is.na(gdpar_mean) && abs(true_val) > 1e-6)
      100 * (gdpar_mean - true_val) / abs(true_val) else NA_real_
    b_dev <- if (!is.na(brms_mean) && abs(true_val) > 1e-6)
      100 * (brms_mean - true_val) / abs(true_val) else NA_real_
    push(sprintf("%-6d %-4s %10.4f %10.4f %10.4f %10.2f %10.2f",
                  k, var_name, true_val, gdpar_mean, brms_mean,
                  g_dev, b_dev))
  }
}
push("```")
push("")

## -- hierarchical scales -------------------------------------------------
push("## gdpar hierarchical scale posteriors")
push("")
stan_fit <- fit_g$fit$fit
stan_vars <- stan_fit$metadata()$stan_variables
for (vn in c("sigma_a", "sigma_b", "sigma_W", "sigma_y")) {
  if (vn %in% stan_vars) {
    push("### ", vn)
    s <- stan_fit$summary(vn)
    push("```")
    capture.output(print(s)) |> (\(z) {
      for (ln in z) push(ln)
    })()
    push("```")
    push("")
  } else {
    push("### ", vn, "  (not present in fit)")
    push("")
  }
}

push("## theta_ref posterior vs true")
push("")
if ("theta_ref" %in% stan_vars) {
  th <- stan_fit$summary("theta_ref")
  push("```")
  capture.output(print(th)) |> (\(z) {
    for (ln in z) push(ln)
  })()
  push("```")
}
push("")
push("True theta_ref: ", paste(true_params$theta_ref, collapse = ", "))
push("")

## -- pointwise ELPD per coord -------------------------------------------
push("## Pointwise ELPD per coordinate")
push("")
push("Decomposes the ELPD gap into per-coordinate contributions. If the")
push("88-pt deficit lives in coord 1 (high-info beta=(1.0,-0.8)), the")
push("hypothesis is hierarchical-prior shrinkage on coord with larger ")
push("posterior magnitudes. If it lives in coords 2-3 (low-info), the")
push("hypothesis flips toward sigma_y miscalibration or pooling pull.")
push("")

draws_arr <- stan_fit$draws(variables = "log_lik", format = "draws_array")
vars <- posterior::variables(draws_arr)
mat <- unclass(posterior::as_draws_matrix(draws_arr))
m <- regexec("^log_lik\\[(\\d+),(\\d+)\\]$", vars)
parsed <- regmatches(vars, m)
ij <- do.call(rbind, lapply(parsed, function(z) as.integer(z[2:3])))
colnames(ij) <- c("i", "k")
n_train <- max(ij[, "i"])

gdpar_elpd_per_k <- numeric(p)
gdpar_se_per_k   <- numeric(p)
for (k in seq_len(p)) {
  cols_k <- which(ij[, "k"] == k)
  ll_k <- mat[, cols_k, drop = FALSE]
  S <- nrow(ll_k); n <- ncol(ll_k)
  log_S <- log(S)
  pointwise <- vapply(seq_len(n), function(i) {
    log_sum_exp(ll_k[, i]) - log_S
  }, numeric(1L))
  gdpar_elpd_per_k[k] <- sum(pointwise)
  gdpar_se_per_k[k]   <- sqrt(n) * stats::sd(pointwise)
}

brms_ll <- brms::log_lik(fit_b$fit)
push("brms log_lik dim: ", paste(dim(brms_ll), collapse = " x "))
brms_elpd_per_k <- numeric(p)
brms_se_per_k   <- numeric(p)
if (length(dim(brms_ll)) == 3L) {
  for (k in seq_len(p)) {
    ll_k <- brms_ll[, , k]
    S <- nrow(ll_k); n <- ncol(ll_k)
    log_S <- log(S)
    pointwise <- vapply(seq_len(n), function(i) {
      log_sum_exp(ll_k[, i]) - log_S
    }, numeric(1L))
    brms_elpd_per_k[k] <- sum(pointwise)
    brms_se_per_k[k]   <- sqrt(n) * stats::sd(pointwise)
  }
} else {
  resp_names <- colnames(brms_ll)
  push("brms log_lik 2D; column-name pattern:")
  push(paste(head(resp_names, 6), collapse = " | "))
  for (k in seq_len(p)) {
    pattern <- sprintf("y%d", k)
    sel <- grep(pattern, resp_names)
    if (length(sel) == 0L) {
      step <- ncol(brms_ll) / p
      sel <- seq((k - 1L) * step + 1L, k * step)
    }
    ll_k <- brms_ll[, sel]
    S <- nrow(ll_k); n <- ncol(ll_k)
    log_S <- log(S)
    pointwise <- vapply(seq_len(n), function(i) {
      log_sum_exp(ll_k[, i]) - log_S
    }, numeric(1L))
    brms_elpd_per_k[k] <- sum(pointwise)
    brms_se_per_k[k]   <- sqrt(n) * stats::sd(pointwise)
  }
}

push("")
push("```")
push(sprintf("%-6s %12s %10s %12s %10s %10s",
              "coord", "gdpar_elpd", "gdpar_se", "brms_elpd",
              "brms_se", "delta"))
for (k in seq_len(p)) {
  delta_k <- gdpar_elpd_per_k[k] - brms_elpd_per_k[k]
  push(sprintf("%-6d %12.2f %10.2f %12.2f %10.2f %10.2f",
                k, gdpar_elpd_per_k[k], gdpar_se_per_k[k],
                brms_elpd_per_k[k], brms_se_per_k[k], delta_k))
}
push(sprintf("%-6s %12.2f %10s %12.2f %10s %10.2f",
              "TOTAL", sum(gdpar_elpd_per_k), "",
              sum(brms_elpd_per_k), "",
              sum(gdpar_elpd_per_k) - sum(brms_elpd_per_k)))
push("```")
push("")

## -- optional wider-prior refit ------------------------------------------
if (identical(Sys.getenv("GDPAR_S7_REFIT_WIDER"), "1")) {
  push("## Refit with wider sigma_b prior")
  push("")
  push("Replacing sigma_b prior with `student_t(3, 0, 5)` (default scale 1).")
  push("If the ELPD gap shrinks, hypothesis confirmed.")
  push("")
  wide_prior <- gdpar_prior(sigma_b = "student_t(3, 0, 5)")
  t2 <- Sys.time()
  iter <- bench_default_iter()
  fit_g_wide <- suppressMessages(suppressWarnings(
    gdpar(
      formula = stats::as.formula(
        sprintf("y ~ %s",
                paste(sc$scenario_meta$gdpar_rhs, collapse = " + "))),
      family = sc$family_spec$gdpar_family,
      amm = sc$scenario_meta$gdpar_amm,
      data = local({
        d <- sc$train
        d$y <- as.matrix(d[, paste0("y", seq_len(p)), drop = FALSE])
        d
      }),
      prior = wide_prior,
      chains = iter$chains,
      iter_warmup = iter$iter_warmup,
      iter_sampling = iter$iter_sampling,
      seed = iter$seed,
      refresh = 0,
      parametrization = "auto",
      skip_id_check = FALSE,
      verbose = FALSE
    )
  ))
  wide_wall <- as.numeric(difftime(Sys.time(), t2, units = "secs"))
  push("gdpar (wide sigma_b) wall_seconds: ", round(wide_wall, 1))
  push("")

  loo_wide <- gdpar_loo(fit_g_wide, aggregation = "subject")
  push("gdpar (wide sigma_b) elpd_loo (subject): ",
       formatC(loo_wide$estimates["elpd_loo", "Estimate"],
                 digits = 2, format = "f"),
       " (se ",
       formatC(loo_wide$estimates["elpd_loo", "SE"],
                 digits = 2, format = "f"), ")")
  push("")
  s_wide <- fit_g_wide$fit$summary(c("sigma_b", "sigma_a"))
  push("```")
  capture.output(print(s_wide)) |> (\(z) {
    for (ln in z) push(ln)
  })()
  push("```")
  push("")
}

## -- write report --------------------------------------------------------
writeLines(report, report_path)
cat("\nReport written to: ", report_path, "\n", sep = "")

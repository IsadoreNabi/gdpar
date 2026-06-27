## calibrate_cp_ncp_multi.R
##
## Multi-coord (p > 1) calibration of preflight_parametrization_multi():
## info-ratio z-test per coord, attribution per coord, e-bfmi global
## (Path B' multi, Phase H.1-H.2 contract). Eight canonical scenarios
## span uniform / heterogeneous / borderline regimes for both a-only
## and W-active families at p = 2 and p = 3.
##
## Usage (from the package root):
##   Rscript inst/benchmarks/calibrate_cp_ncp_multi.R
##
## The script loads the package via devtools::load_all() and writes
## inst/benchmarks/results/cp_ncp_hit_rate_multi.csv with the long-tidy
## per-(scenario, k, component) table. It prints a summary hit rate
## table at the end.
##
## Wall time: approximately 15-25 minutes on the reference machine
## (1 Stan compile multi + 8 preflights of 200/200 each + 8
## contrastive short fits with the alternative parametrization).
##
## Scenarios
##  m1  p=2  uniform high-info, a-only          expected (TRUE , TRUE )
##  m2  p=2  uniform low-info , a-only          expected (FALSE, FALSE)
##  m3  p=2  heterogeneous mixto, a-only        expected (TRUE , FALSE)
##  m4  p=2  W high-info alone, no a            expected cp_W = TRUE
##  m5  p=2  a high + W high (confounded)       expected a=(TRUE,TRUE)
##  m6  p=2  borderline tau_cp zone, a-only     expected NA (post-hoc)
##  m7  p=3  uniform high-info, a-only          expected (TRUE,TRUE,TRUE)
##  m8  p=3  heterogeneous a-only               expected (TRUE,FALSE,TRUE)

suppressPackageStartupMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  library(devtools)
  load_all(".", quiet = TRUE)
})

# --- helpers ----------------------------------------------------------

simulate_a_block <- function(n, p, theta_ref_true, beta_mat,
                             x_mat, sigma_vec) {
  stopifnot(nrow(beta_mat) == p, ncol(beta_mat) == ncol(x_mat))
  eta <- matrix(theta_ref_true, n, p, byrow = TRUE) +
    x_mat %*% t(beta_mat)
  if (length(sigma_vec) == 1L) {
    sigma_vec <- rep(sigma_vec, p)
  }
  noise <- vapply(seq_len(p), function(k) rnorm(n, 0, sigma_vec[k]),
                  numeric(n))
  eta + noise
}

simulate_W_block_per_k <- function(n, p, theta_ref_true,
                                    w_deg1_per_k, w_deg2_per_k,
                                    x_mat) {
  eta_W <- matrix(0, n, p)
  for (k in seq_len(p)) {
    c1 <- theta_ref_true[k]     * w_deg1_per_k[[k]]
    c2 <- theta_ref_true[k]^2   * w_deg2_per_k[[k]]
    eta_W[, k] <- drop(x_mat %*% (c1 + c2))
  }
  eta_W
}

prepare_stan_data <- function(spec, df, formula_rhs, y_mat,
                              theta_anchor = 0) {
  fam <- gdpar_family_multi("gaussian", p = spec$p)
  if (!is.null(spec$W)) {
    spec$W <- materialize_W_basis(spec$W, p = spec$p)
  }
  design <- build_amm_design(spec, df, formula_rhs = formula_rhs)
  stan_data <- assemble_stan_data(design, fam, spec, y_mat,
                                  theta_anchor = theta_anchor)
  list(spec = spec, design = design, stan_data = stan_data, family = fam)
}

run_preflight <- function(prep, prior, seed = 42L) {
  preflight_parametrization_multi(
    prior, prep$stan_data, prep$spec,
    preflight_seed = seed, verbose = FALSE,
    aggregation = "any_ncp"
  )
}

run_contrastive_fit <- function(prep, prior, cp_a, cp_W,
                                cp_a_per_k = NULL, seed = 42L) {
  stan_src <- generate_stan_code_multi(
    prior, cp_a = isTRUE(cp_a), cp_W = isTRUE(cp_W),
    cp_a_per_k = cp_a_per_k
  )
  stan_path <- write_stan_to_tempfile(stan_src)
  cs_model <- cmdstanr::cmdstan_model(stan_path)
  sd <- prep$stan_data
  if (!is.null(cp_a_per_k) && length(unique(cp_a_per_k)) > 1L) {
    sd$cp_a_per_k_data <- as.integer(cp_a_per_k)
  }
  fit <- cs_model$sample(
    data = sd, chains = 2L, parallel_chains = 1L,
    iter_warmup = 200L, iter_sampling = 200L,
    refresh = 0L, seed = seed,
    show_messages = FALSE, show_exceptions = FALSE
  )
  diag <- fit$diagnostic_summary(quiet = TRUE)
  list(
    n_divergent = as.integer(sum(diag$num_divergent)),
    ebfmi_min   = min(diag$ebfmi),
    div_pct     = sum(diag$num_divergent) / (length(diag$num_divergent) * 200L)
  )
}

# --- scenario constructors -------------------------------------------

scenarios <- list()

scenarios[["m1_homog_high_p2"]] <- list(
  description = "p=2, both coords high-info, a-only (no W).",
  p = 2L,
  expected_cp_a_per_k = c(TRUE, TRUE),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 500L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_mat <- cbind(x1, x2, x3, x4)
    beta_mat <- rbind(c( 1.0, -0.8,  0.6, -0.4),
                       c( 0.9, -0.7,  0.5, -0.3))
    y_mat <- simulate_a_block(n, 2L, c(0.5, -0.5), beta_mat, x_mat,
                              sigma_vec = c(0.10, 0.10))
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 2L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4)))
  }
)

scenarios[["m2_homog_low_p2"]] <- list(
  description = "p=2, both coords low-info, a-only (no W).",
  p = 2L,
  expected_cp_a_per_k = c(FALSE, FALSE),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 40L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_mat <- cbind(x1, x2, x3, x4)
    beta_mat <- rbind(c( 0.02, -0.02,  0.01, -0.01),
                       c( 0.02, -0.02,  0.01, -0.01))
    y_mat <- simulate_a_block(n, 2L, c(0.5, -0.5), beta_mat, x_mat,
                              sigma_vec = c(3.0, 3.0))
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 2L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4)))
  }
)

scenarios[["m3_heterog_mixto_p2"]] <- list(
  description = "p=2, k=1 high-info, k=2 low-info, a-only (no W).",
  p = 2L,
  expected_cp_a_per_k = c(TRUE, FALSE),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 200L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_mat <- cbind(x1, x2, x3, x4)
    beta_mat <- rbind(c( 1.0, -0.8,  0.6, -0.4),
                       c( 0.02, -0.02,  0.01, -0.01))
    y_mat <- simulate_a_block(n, 2L, c(0.5, -0.5), beta_mat, x_mat,
                              sigma_vec = c(0.10, 2.0))
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 2L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4)))
  }
)

scenarios[["m4_W_high_only_p2"]] <- list(
  description = "p=2, W alone, both coords high-info via W (no a).",
  p = 2L,
  expected_cp_a_per_k = c(NA, NA),
  expected_cp_W       = TRUE,
  expected_cp_W_per_k = c(TRUE, TRUE),
  generator = function(seed) {
    set.seed(seed)
    n <- 500L
    x1 <- rnorm(n); x2 <- rnorm(n)
    x_mat <- cbind(x1, x2)
    eta_W <- simulate_W_block_per_k(
      n, 2L, c(0.3, -0.3),
      w_deg1_per_k = list(c(2.0, -1.5), c(1.8, -1.4)),
      w_deg2_per_k = list(c(-1.0, 0.8), c(-0.9, 0.7)),
      x_mat = x_mat
    )
    y_mat <- matrix(c(0.3, -0.3), n, 2L, byrow = TRUE) + eta_W +
      matrix(rnorm(n * 2L, 0, 0.05), n, 2L)
    df <- data.frame(x1 = x1, x2 = x2)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2,
         spec = amm_spec(p = 2L,
                          dims = dimwise(a = NULL, b = NULL),
                          W = W_basis(type = "polynomial", degree = 2),
                          x_vars = c("x1", "x2")))
  }
)

scenarios[["m5_both_high_W_active_p2"]] <- list(
  description = "p=2, a high + W high (structural confounding via overlapping x_vars).",
  p = 2L,
  expected_cp_a_per_k = c(TRUE, TRUE),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 500L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_a <- cbind(x1, x2, x3, x4)
    x_W <- cbind(x1, x2)
    beta_mat <- rbind(c( 1.0, -0.8,  0.6, -0.4),
                       c( 0.9, -0.7,  0.5, -0.3))
    eta_a <- x_a %*% t(beta_mat)
    eta_W <- simulate_W_block_per_k(
      n, 2L, c(0.3, -0.3),
      w_deg1_per_k = list(c(2.0, -1.5), c(1.8, -1.4)),
      w_deg2_per_k = list(c(-1.0, 0.8), c(-0.9, 0.7)),
      x_mat = x_W
    )
    y_mat <- matrix(c(0.3, -0.3), n, 2L, byrow = TRUE) +
      eta_a + eta_W + matrix(rnorm(n * 2L, 0, 0.10), n, 2L)
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 2L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4),
                          W = W_basis(type = "polynomial", degree = 2),
                          x_vars = c("x1", "x2")))
  }
)

scenarios[["m6_borderline_p2"]] <- list(
  description = "p=2, borderline info-ratio zone (~2-5) for a-only.",
  p = 2L,
  expected_cp_a_per_k = c(NA, NA),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 100L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_mat <- cbind(x1, x2, x3, x4)
    beta_mat <- rbind(c( 0.30, -0.25,  0.20, -0.15),
                       c( 0.28, -0.23,  0.18, -0.13))
    y_mat <- simulate_a_block(n, 2L, c(0.5, -0.5), beta_mat, x_mat,
                              sigma_vec = c(0.5, 0.5))
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 2L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4)))
  }
)

scenarios[["m7_homog_high_p3"]] <- list(
  description = "p=3, uniform high-info, a-only (no W).",
  p = 3L,
  expected_cp_a_per_k = c(TRUE, TRUE, TRUE),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 500L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_mat <- cbind(x1, x2, x3, x4)
    beta_mat <- rbind(c( 1.0, -0.8,  0.6, -0.4),
                       c( 0.9, -0.7,  0.5, -0.3),
                       c( 0.8, -0.6,  0.4, -0.2))
    y_mat <- simulate_a_block(n, 3L, c(0.5, -0.5, 0.0), beta_mat,
                              x_mat, sigma_vec = c(0.10, 0.10, 0.10))
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 3L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4)))
  }
)

scenarios[["m8_heterog_p3"]] <- list(
  description = "p=3, k=1 high, k=2 low, k=3 high, a-only.",
  p = 3L,
  expected_cp_a_per_k = c(TRUE, FALSE, TRUE),
  expected_cp_W       = NA,
  expected_cp_W_per_k = c(NA, NA, NA),
  generator = function(seed) {
    set.seed(seed)
    n <- 300L
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
    x_mat <- cbind(x1, x2, x3, x4)
    beta_mat <- rbind(c( 1.0, -0.8,  0.6, -0.4),
                       c( 0.02, -0.02,  0.01, -0.01),
                       c( 0.9, -0.7,  0.5, -0.3))
    y_mat <- simulate_a_block(n, 3L, c(0.5, -0.5, 0.0), beta_mat,
                              x_mat, sigma_vec = c(0.10, 2.0, 0.10))
    df <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    list(df = df, y_mat = y_mat,
         formula_rhs = ~ x1 + x2 + x3 + x4,
         spec = amm_spec(p = 3L,
                          dims = dimwise(a = ~ x1 + x2 + x3 + x4)))
  }
)

# --- main loop --------------------------------------------------------

prior <- gdpar_prior()
rows <- list()

for (s_name in names(scenarios)) {
  s <- scenarios[[s_name]]
  cat(sprintf("\n[%s] %s\n", s_name, s$description))
  d <- s$generator(seed = 42L)
  prep <- prepare_stan_data(d$spec, d$df, d$formula_rhs, d$y_mat,
                            theta_anchor = 0)
  cat("  running preflight (multi NCP short)...\n")
  pre <- run_preflight(prep, prior, seed = 42L)
  per_dim <- pre$report$per_dim
  cp_a_pred_per_k <- pre$cp_a_per_k
  cp_W_pred_per_k <- pre$cp_W_per_k
  n_div_pred <- per_dim$n_divergent[1]
  ebfmi_min  <- per_dim$ebfmi_min[1]
  cat(sprintf("  preflight done: cp_a_per_k=(%s) cp_W_per_k=(%s) ndiv=%d ebfmi_min=%.3f\n",
              paste(cp_a_pred_per_k, collapse = ","),
              paste(cp_W_pred_per_k, collapse = ","),
              as.integer(n_div_pred), as.numeric(ebfmi_min)))

  alt_cp_a_per_k <- !cp_a_pred_per_k
  alt_cp_W <- !isTRUE(pre$cp_W)
  if (length(unique(alt_cp_a_per_k)) > 1L) {
    alt_cp_a_global <- TRUE
  } else {
    alt_cp_a_global <- alt_cp_a_per_k[1]
  }
  cat("  running contrastive fit (alternative parametrization)...\n")
  alt <- run_contrastive_fit(
    prep, prior,
    cp_a = alt_cp_a_global, cp_W = alt_cp_W,
    cp_a_per_k = if (length(unique(alt_cp_a_per_k)) > 1L) alt_cp_a_per_k else NULL,
    seed = 42L
  )
  cat(sprintf("  contrastive done: ndiv_alt=%d ebfmi_alt=%.3f\n",
              alt$n_divergent, alt$ebfmi_min))

  for (k in seq_len(s$p)) {
    row_a <- per_dim[per_dim$component == "a" & per_dim$dim == k, ]
    truth_a <- s$expected_cp_a_per_k[k]
    pred_a <- cp_a_pred_per_k[k]
    hit_a <- if (is.na(truth_a) || is.na(pred_a)) NA
             else identical(truth_a, pred_a)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario        = s_name,
      p               = s$p,
      k               = k,
      component       = "a",
      regime_truth    = truth_a,
      regime_pred     = pred_a,
      hit             = hit_a,
      decision_reason = if (nrow(row_a) == 1L) row_a$decision_reason
                        else NA_character_,
      n_div_pred      = n_div_pred,
      n_div_alt       = alt$n_divergent,
      ebfmi_min       = round(ebfmi_min, 3),
      t_attr          = if (nrow(row_a) == 1L) round(row_a$t_attribution, 3) else NA_real_,
      t_info_cp       = if (nrow(row_a) == 1L) round(row_a$t_info_cp, 3) else NA_real_,
      t_info_ncp      = if (nrow(row_a) == 1L) round(row_a$t_info_ncp, 3) else NA_real_,
      stringsAsFactors = FALSE
    )
    row_W <- per_dim[per_dim$component == "W" & per_dim$dim == k, ]
    truth_W <- s$expected_cp_W_per_k[k]
    pred_W <- cp_W_pred_per_k[k]
    hit_W <- if (is.na(truth_W) || is.na(pred_W)) NA
             else identical(truth_W, pred_W)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario        = s_name,
      p               = s$p,
      k               = k,
      component       = "W",
      regime_truth    = truth_W,
      regime_pred     = pred_W,
      hit             = hit_W,
      decision_reason = if (nrow(row_W) == 1L) row_W$decision_reason
                        else NA_character_,
      n_div_pred      = n_div_pred,
      n_div_alt       = alt$n_divergent,
      ebfmi_min       = round(ebfmi_min, 3),
      t_attr          = if (nrow(row_W) == 1L) round(row_W$t_attribution, 3) else NA_real_,
      t_info_cp       = if (nrow(row_W) == 1L) round(row_W$t_info_cp, 3) else NA_real_,
      t_info_ncp      = if (nrow(row_W) == 1L) round(row_W$t_info_ncp, 3) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, rows)

out_dir <- "inst/benchmarks/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, "cp_ncp_hit_rate_multi.csv")
write.csv(results, out_csv, row.names = FALSE)

cat("\n\n=== Full results (long-tidy) ===\n")
print(results, row.names = FALSE, max = 1000)

cat("\n=== Hit rate per component (truths only) ===\n")
sub <- results[!is.na(results$regime_truth) & !is.na(results$hit), ]
if (nrow(sub) > 0L) {
  hits <- aggregate(hit ~ component, data = sub,
                    FUN = function(x) mean(as.logical(x), na.rm = TRUE))
  print(hits, row.names = FALSE)
} else {
  cat("(no rows with truth defined)\n")
}

cat(sprintf("\nWrote: %s\n", out_csv))
cat("Done.\n")

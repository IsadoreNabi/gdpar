## Synthetic adversarial scenarios for the multivariate validation
## suite (Block 6, coord-wise factorization).
##
## Nine generators S0..S8, each returns a list with:
##   $train, $test:        data.frame splits (n_train + n_test = 1000)
##   $true_params:         named list of ground-truth parameters
##   $family_spec:         per-method family objects (gdpar, brms, mgcv,
##                          INLA, rstanarm)
##   $scenario_meta:       per-method formulas, gdpar amm_spec object,
##                          tier (piso / sanity / adversarial / stress),
##                          and the scenario_id string used in the CSV.
##
## All generators use a deterministic seed = 42L + scenario_offset to
## make individual runs reproducible.
##
## Family policy: gdpar_family_multi() currently rejects heterogeneous
## families across coordinates ("gdpar_unsupported_feature_error",
## families.R lines 438-449). Until Phase F.x lifts that restriction,
## S5 is implemented as homogeneous Bernoulli p=2 with one high-info
## and one low-info coordinate, which exercises the classification
## branch (AUC / log_loss / Brier) of the bench harness without
## violating the API.
##
## n_train = 800, n_test = 200 per scenario; total 1000 rows.

scenarios_synthetic_hard_multi <- function() {
  list(
    S0 = gen_S0_calibrate_replica,
    S1 = gen_S1_sanity_disjoint,
    S2 = gen_S2_W_only,
    S3 = gen_S3_overlap_controlled,
    S4 = gen_S4_near_collinear,
    S5 = gen_S5_bernoulli_heterog,
    S6 = gen_S6_C4bis_overlap_stress,
    S7 = gen_S7_heterog_p3,
    S8 = gen_S8_residual_correlation
  )
}

.split_train_test <- function(full, n_train, n_test) {
  list(train = full[seq_len(n_train), , drop = FALSE],
       test  = full[(n_train + 1L):(n_train + n_test), , drop = FALSE])
}

.family_spec_gaussian_homog <- function(p) {
  list(
    gdpar_family    = gdpar_family_multi("gaussian", p = p),
    brms_family     = "gaussian",
    mgcv_family     = replicate(p, stats::gaussian(), simplify = FALSE),
    inla_family     = rep("gaussian", p),
    rstanarm_family = replicate(p, stats::gaussian(), simplify = FALSE)
  )
}

.family_spec_bernoulli_homog <- function(p) {
  list(
    gdpar_family    = gdpar_family_multi("bernoulli", p = p),
    brms_family     = "bernoulli",
    mgcv_family     = replicate(p, stats::binomial(link = "logit"),
                                 simplify = FALSE),
    inla_family     = rep("binomial", p),
    rstanarm_family = replicate(p, stats::binomial(link = "logit"),
                                 simplify = FALSE)
  )
}

.bench_formula_per_coord <- function(rhs_terms, p) {
  rhs <- paste(rhs_terms, collapse = " + ")
  lapply(seq_len(p), function(k) {
    stats::as.formula(sprintf("y%d ~ %s", k, rhs))
  })
}

gen_S0_calibrate_replica <- function(seed = 42L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.5, -0.5)
  beta_a_true <- matrix(c(0.5, -0.5, -0.4, 0.4), nrow = 2L, ncol = p)
  sigma_y_true <- c(0.3, 0.3)
  eta <- cbind(
    theta_ref_true[1] + beta_a_true[1, 1] * x1 + beta_a_true[2, 1] * x2,
    theta_ref_true[2] + beta_a_true[1, 2] * x1 + beta_a_true[2, 2] * x2
  )
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S0_calibrate_replica",
      tier        = "piso",
      p           = p,
      gdpar_amm   = amm_spec(p = p,
                              dims = dimwise(a = ~ x1 + x2, b = NULL)),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2")
    )
  )
}

gen_S1_sanity_disjoint <- function(seed = 43L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n)
  theta_ref_true <- c(0.3, -0.4)
  beta_a_true <- matrix(c(0.6, -0.4, -0.5, 0.5), nrow = 2L, ncol = p)
  sigma_y_true <- c(0.4, 0.4)
  eta <- cbind(
    theta_ref_true[1] + beta_a_true[1, 1] * x1 + beta_a_true[2, 1] * x2,
    theta_ref_true[2] + beta_a_true[1, 2] * x1 + beta_a_true[2, 2] * x2
  )
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L],
                      x1 = x1, x2 = x2, z1 = z1, z2 = z2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S1_sanity_disjoint",
      tier        = "sanity",
      p           = p,
      gdpar_amm   = amm_spec(
        p = p,
        dims = dimwise(a = ~ x1 + x2, b = NULL),
        W = W_basis(type = "polynomial", degree = 2L),
        x_vars = c("z1", "z2")
      ),
      mgcv_formula_per_coord = .bench_formula_per_coord(
        c("x1", "x2", "z1", "z2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(
        c("x1", "x2", "z1", "z2"), p),
      rstanarm_formula_per_coord = .bench_formula_per_coord(
        c("x1", "x2", "z1", "z2"), p),
      brms_rhs    = "x1 + x2 + z1 + z2",
      gdpar_rhs   = c("x1", "x2", "z1", "z2")
    )
  )
}

gen_S2_W_only <- function(seed = 44L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  z1 <- rnorm(n); z2 <- rnorm(n)
  theta_ref_true <- c(0.4, -0.3)
  W_true <- matrix(c(0.5, -0.3, -0.4, 0.4), nrow = 2L, ncol = 2L)
  sigma_y_true <- c(0.4, 0.4)
  eta <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(p)) {
    eta[, k] <- theta_ref_true[k] +
      (theta_ref_true[k]^2 - 0^2) *
        (W_true[1L, 1L] * z1 + W_true[1L, 2L] * z2) +
      (theta_ref_true[k]   - 0)   *
        (W_true[2L, 1L] * z1 + W_true[2L, 2L] * z2)
  }
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], z1 = z1, z2 = z2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        W = W_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S2_W_only",
      tier        = "adversarial",
      p           = p,
      gdpar_amm   = amm_spec(
        p = p,
        dims = dimwise(a = NULL, b = NULL),
        W = W_basis(type = "polynomial", degree = 2L),
        x_vars = c("z1", "z2")
      ),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("z1", "z2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("z1", "z2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("z1", "z2"), p),
      brms_rhs    = "z1 + z2",
      gdpar_rhs   = c("z1", "z2")
    )
  )
}

gen_S3_overlap_controlled <- function(seed = 45L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.5, -0.5)
  beta_a_true <- matrix(c(0.5, -0.3, -0.4, 0.4), nrow = 2L, ncol = p)
  W_true <- matrix(c(0.2, -0.1, -0.15, 0.18), nrow = 2L, ncol = 2L)
  sigma_y_true <- c(0.35, 0.35)
  eta <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(p)) {
    eta[, k] <- theta_ref_true[k] +
      beta_a_true[1L, k] * x1 + beta_a_true[2L, k] * x2 +
      (theta_ref_true[k]^2) *
        (W_true[1L, 1L] * x1 + W_true[1L, 2L] * x2) +
      (theta_ref_true[k]) *
        (W_true[2L, 1L] * x1 + W_true[2L, 2L] * x2)
  }
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        W = W_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S3_overlap_controlled",
      tier        = "adversarial",
      p           = p,
      rigor       = "fast",
      gdpar_amm   = amm_spec(
        p = p,
        dims = dimwise(a = ~ x1 + x2, b = NULL),
        W = W_basis(type = "polynomial", degree = 2L),
        x_vars = c("x1", "x2")
      ),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2")
    )
  )
}

gen_S4_near_collinear <- function(seed = 46L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.05)
  theta_ref_true <- c(0.4, -0.4)
  beta_a_true <- matrix(c(0.6, -0.4, -0.5, 0.5), nrow = 2L, ncol = p)
  sigma_y_true <- c(0.5, 0.5)
  eta <- cbind(
    theta_ref_true[1] + beta_a_true[1, 1] * x1 + beta_a_true[2, 1] * x2,
    theta_ref_true[2] + beta_a_true[1, 2] * x1 + beta_a_true[2, 2] * x2
  )
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S4_near_collinear",
      tier        = "adversarial",
      p           = p,
      gdpar_amm   = amm_spec(p = p,
                              dims = dimwise(a = ~ x1 + x2, b = NULL)),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2")
    )
  )
}

gen_S5_bernoulli_heterog <- function(seed = 47L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.0, 0.0)
  beta_a_high <- c(1.5, -1.2)
  beta_a_low  <- c(0.3, -0.2)
  eta <- cbind(
    theta_ref_true[1] + beta_a_high[1] * x1 + beta_a_high[2] * x2,
    theta_ref_true[2] + beta_a_low[1]  * x1 + beta_a_low[2]  * x2
  )
  prob <- 1 / (1 + exp(-eta))
  Y <- matrix(stats::rbinom(n * p, size = 1L, prob = as.vector(prob)),
              nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = list(high = beta_a_high,
                                       low = beta_a_low)),
    family_spec = .family_spec_bernoulli_homog(p),
    scenario_meta = list(
      scenario_id = "S5_bernoulli_heterog",
      tier        = "adversarial",
      p           = p,
      gdpar_amm   = amm_spec(p = p,
                              dims = dimwise(a = ~ x1 + x2, b = NULL)),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2"),
      classification_metrics = TRUE
    )
  )
}

gen_S6_C4bis_overlap_stress <- function(seed = 48L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.6, -0.6)
  beta_a_true <- matrix(c(0.8, -0.5, -0.7, 0.6), nrow = 2L, ncol = p)
  W_true <- matrix(c(0.3, -0.2, -0.2, 0.3), nrow = 2L, ncol = 2L)
  sigma_y_true <- c(0.4, 0.4)
  eta <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(p)) {
    eta[, k] <- theta_ref_true[k] +
      beta_a_true[1L, k] * x1 + beta_a_true[2L, k] * x2 +
      (theta_ref_true[k]^2) *
        (W_true[1L, 1L] * x1 + W_true[1L, 2L] * x2) +
      (theta_ref_true[k]) *
        (W_true[2L, 1L] * x1 + W_true[2L, 2L] * x2)
  }
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        W = W_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S6_C4bis_overlap_stress",
      tier        = "stress",
      p           = p,
      rigor       = "fast",
      gdpar_amm   = amm_spec(
        p = p,
        dims = dimwise(a = ~ x1 + x2, b = NULL),
        W = W_basis(type = "polynomial", degree = 2L),
        x_vars = c("x1", "x2")
      ),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2")
    )
  )
}

gen_S7_heterog_p3 <- function(seed = 49L) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 3L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.5, -0.3, 0.2)
  beta_a_true <- matrix(c(1.0, -0.8,
                            0.2, -0.15,
                            0.18, -0.12),
                         nrow = 2L, ncol = p)
  sigma_y_true <- c(0.3, 0.5, 0.5)
  eta <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(p)) {
    eta[, k] <- theta_ref_true[k] +
      beta_a_true[1L, k] * x1 + beta_a_true[2L, k] * x2
  }
  Y <- eta + matrix(rnorm(n * p, sd = rep(sigma_y_true, each = n)),
                    nrow = n, ncol = p)
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], y3 = Y[, 3L],
                      x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        sigma_y = sigma_y_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S7_heterog_p3",
      tier        = "stress",
      p           = p,
      gdpar_amm   = amm_spec(p = p,
                              dims = dimwise(a = ~ x1 + x2, b = NULL)),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2")
    )
  )
}

gen_S8_residual_correlation <- function(seed = 50L,
                                          rho_true = 0.5,
                                          sigma_y = 0.4) {
  set.seed(seed)
  n_train <- 800L; n_test <- 200L; n <- n_train + n_test; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.4, -0.4)
  beta_a_true <- matrix(c(0.6, -0.4, -0.5, 0.5), nrow = 2L, ncol = p)
  eta <- cbind(
    theta_ref_true[1] + beta_a_true[1, 1] * x1 + beta_a_true[2, 1] * x2,
    theta_ref_true[2] + beta_a_true[1, 2] * x1 + beta_a_true[2, 2] * x2
  )
  Sigma <- matrix(c(1, rho_true, rho_true, 1) * sigma_y^2,
                   nrow = 2L, ncol = 2L)
  if (!requireNamespace("mvnfast", quietly = TRUE)) {
    L <- t(chol(Sigma))
    eps_iid <- matrix(rnorm(n * p), nrow = n, ncol = p)
    eps <- eps_iid %*% t(L)
  } else {
    eps <- mvnfast::rmvn(n = n, mu = c(0, 0), sigma = Sigma)
  }
  Y <- eta + eps
  full <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x1 = x1, x2 = x2)
  split <- .split_train_test(full, n_train, n_test)
  list(
    train = split$train, test = split$test,
    true_params = list(theta_ref = theta_ref_true,
                        beta_a = beta_a_true,
                        sigma_y = c(sigma_y, sigma_y),
                        rho_true = rho_true),
    family_spec = .family_spec_gaussian_homog(p),
    scenario_meta = list(
      scenario_id = "S8_residual_correlation",
      tier        = "stress",
      p           = p,
      brms_rescor_recovery = TRUE,
      gdpar_amm   = amm_spec(p = p,
                              dims = dimwise(a = ~ x1 + x2, b = NULL)),
      mgcv_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      inla_formula_per_coord = .bench_formula_per_coord(c("x1", "x2"), p),
      rstanarm_formula_per_coord =
        .bench_formula_per_coord(c("x1", "x2"), p),
      brms_rhs    = "x1 + x2",
      gdpar_rhs   = c("x1", "x2")
    )
  )
}

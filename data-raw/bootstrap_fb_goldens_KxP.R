# Bootstrap script for Sub-bloque 9.3.d (Bloque 9, Sesion B9.5,
# 2026-05-27; deuda D69 substantively closed in Sesion B9.6,
# 2026-05-27) under decisions I.iv + J.iv.A + M.iv + N.a. Parallel to
# data-raw/bootstrap_eb_goldens_D.R (8.6.D EB Path C) but writes FB
# Path C golden snapshots for the K > 1 + p > 1 regime via the new
# piece inst/stan/_canonical_pieces/amm_canonical_pmulti_KxP.stan plus
# the unified assembler .assemble_stan_data_KxP(path = "FB") and the
# driver .gdpar_fb_KxP_fit().
#
# Usage:
#
#   GDPAR_BOOTSTRAP_FB_GOLDENS=1 \
#   GDPAR_BOOTSTRAP_FB_GOLDENS_KxP=1 \
#     Rscript data-raw/bootstrap_fb_goldens_KxP.R
#
# Both env vars must be set: the outer flag GDPAR_BOOTSTRAP_FB_GOLDENS
# scopes any FB bootstrap (forward-compatible with future FB roster
# scripts), and the inner flag GDPAR_BOOTSTRAP_FB_GOLDENS_KxP scopes
# the run to the 9.3.d KxP roster. Without both flags the script is
# a no-op informational message.
#
# Reserved canonical seeds for the K.c minimum roster (DESIGN_9_3_D
# §4 of B9.4 + K.c canonized in B9.5):
#   g_KxP_001  gaussian   K=2  p=2   seed = 91001L
#   g_KxP_002  beta       K=2  p=2   seed = 91002L
#   g_KxP_003  gamma      K=2  p=3   seed = 91003L
#   g_KxP_004  student_t  K=3  p=2   seed = 91004L
# The seeds in 91005..91099 stay reserved for the expanded roster
# (heterogeneous slot + K=3 p=3) deferred to B9.7+ if the bench
# external 9.2 requires it.

if (!nzchar(Sys.getenv("GDPAR_BOOTSTRAP_FB_GOLDENS")) ||
    !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_FB_GOLDENS_KxP"))) {
  message(
    "GDPAR_BOOTSTRAP_FB_GOLDENS / GDPAR_BOOTSTRAP_FB_GOLDENS_KxP not set; ",
    "bootstrap_fb_goldens_KxP is a no-op."
  )
  quit(save = "no", status = 0L)
}

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

stopifnot(requireNamespace("cmdstanr", quietly = TRUE))
stopifnot(requireNamespace("posterior", quietly = TRUE))

set.seed(91001L)
data_dir <- file.path("tests", "testthat", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(data_dir, "golden_manifest.csv")
gdpar:::.gdpar_golden_manifest_init(manifest_path)

R_version_str <- paste0(R.version$major, ".", R.version$minor)
cmdstan_version_str <- gdpar:::.gdpar_cmdstan_version()
bootstrap_date <- format(Sys.Date(), "%Y-%m-%d")

ITER_WARMUP <- 200L
ITER_SAMPLING <- 200L
CHAINS <- 2L

# ---- Path C FB K x p fit closures (Sesion B9.6 D69 substantive closure) ----
#
# Each closure builds a deterministic synthetic dataset matching the
# canonical (K, p, family) tuple of the roster entry and calls the
# unified driver .gdpar_fb_KxP_fit() that composes design + assembler
# (path = "FB") + cmdstanr sampling on the canonical piece. The
# synthetic data is small (n = 80-120) so the bootstrap completes in
# minutes; the goal is the bit-exact comparator template golden, not
# inferential validity.

.fit_fb_KxP_gaussian_K2_p2 <- function(entry) {
  set.seed(entry$seed)
  n <- 80L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 2L)
  for (j in 1:2) {
    mu_j <- c(1.0, -0.5)[j] + c(0.5, -0.3)[j] * d$x1
    ls_j <- c(-0.2, 0.0)[j] + c(0.2, 0.1)[j] * d$x2
    y_mat[, j] <- rnorm(n, mean = mu_j, sd = exp(ls_j))
  }
  amm_list <- list(
    mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
    sigma = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
  )
  gdpar:::.gdpar_fb_KxP_fit(
    y_matrix      = y_mat,
    family        = gdpar_family("gaussian"),
    amm_list      = amm_list,
    data          = d,
    formula_rhs   = ~ x1 + x2,
    iter_warmup   = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    chains        = CHAINS,
    seed          = entry$seed
  )
}

.fit_fb_KxP_beta_K2_p2 <- function(entry) {
  set.seed(entry$seed)
  n <- 80L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 2L)
  for (j in 1:2) {
    mu_eta <- c(0.0, 0.5)[j] + c(0.4, -0.3)[j] * d$x1
    phi_eta <- c(1.5, 1.2)[j] + c(0.2, 0.1)[j] * d$x2
    mu_j <- plogis(mu_eta)
    phi_j <- exp(phi_eta)
    a_par <- mu_j * phi_j
    b_par <- (1 - mu_j) * phi_j
    y_mat[, j] <- pmin(pmax(rbeta(n, a_par, b_par), 1e-3), 1 - 1e-3)
  }
  amm_list <- list(
    mu  = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
    phi = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
  )
  gdpar:::.gdpar_fb_KxP_fit(
    y_matrix      = y_mat,
    family        = gdpar_family("beta"),
    amm_list      = amm_list,
    data          = d,
    formula_rhs   = ~ x1 + x2,
    iter_warmup   = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    chains        = CHAINS,
    seed          = entry$seed
  )
}

.fit_fb_KxP_gamma_K2_p3 <- function(entry) {
  set.seed(entry$seed)
  n <- 100L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 3L)
  for (j in 1:3) {
    mu_eta    <- c(0.5, 1.0, 0.0)[j] + c(0.2, -0.1, 0.1)[j] * d$x1
    shape_eta <- c(1.0, 0.5, 0.3)[j] + c(0.1, 0.0, -0.05)[j] * d$x2
    mu_j    <- exp(mu_eta)
    shape_j <- exp(shape_eta)
    rate_j  <- shape_j / mu_j
    y_mat[, j] <- rgamma(n, shape = shape_j, rate = rate_j)
  }
  amm_list <- list(
    mu    = amm_spec(p = 3L, dims = dimwise(a = ~ x1)),
    shape = amm_spec(p = 3L, dims = dimwise(a = ~ x2))
  )
  gdpar:::.gdpar_fb_KxP_fit(
    y_matrix      = y_mat,
    family        = gdpar_family("gamma"),
    amm_list      = amm_list,
    data          = d,
    formula_rhs   = ~ x1 + x2,
    iter_warmup   = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    chains        = CHAINS,
    seed          = entry$seed
  )
}

.fit_fb_KxP_student_t_K3_p2 <- function(entry) {
  set.seed(entry$seed)
  n <- 80L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 2L)
  for (j in 1:2) {
    mu_j    <- c(1.0, -0.5)[j] + c(0.3, -0.2)[j] * d$x1
    sigma_j <- exp(c(-0.5, -0.3)[j] + c(0.1, 0.05)[j] * d$x2)
    nu_j    <- exp(c(1.5, 1.8)[j])
    y_mat[, j] <- mu_j + sigma_j * rt(n, df = nu_j)
  }
  amm_list <- list(
    mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
    sigma = amm_spec(p = 2L, dims = dimwise(a = ~ x2)),
    nu    = amm_spec(p = 2L, dims = dimwise(a = ~ 1))
  )
  gdpar:::.gdpar_fb_KxP_fit(
    y_matrix      = y_mat,
    family        = gdpar_family("student_t"),
    amm_list      = amm_list,
    data          = d,
    formula_rhs   = ~ x1 + x2,
    iter_warmup   = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    chains        = CHAINS,
    seed          = entry$seed
  )
}

.fb_KxP_roster <- list(
  list(name = "g_KxP_001", family = "gaussian",  K = 2L, p = 2L,
       seed = 91001L, basis_type = "polynomial",
       fit_fn = .fit_fb_KxP_gaussian_K2_p2),
  list(name = "g_KxP_002", family = "beta",      K = 2L, p = 2L,
       seed = 91002L, basis_type = "polynomial",
       fit_fn = .fit_fb_KxP_beta_K2_p2),
  list(name = "g_KxP_003", family = "gamma",     K = 2L, p = 3L,
       seed = 91003L, basis_type = "polynomial",
       fit_fn = .fit_fb_KxP_gamma_K2_p3),
  list(name = "g_KxP_004", family = "student_t", K = 3L, p = 2L,
       seed = 91004L, basis_type = "polynomial",
       fit_fn = .fit_fb_KxP_student_t_K3_p2)
)

manifest_rows_written <- 0L

for (entry in .fb_KxP_roster) {
  rds_name <- sprintf("golden_fb_KxP_%s_%s_K%d_p%d.rds",
                       entry$family, entry$basis_type, entry$K, entry$p)
  rds_path <- file.path(data_dir, rds_name)
  rds_path_manifest <- file.path("data", rds_name)
  if (file.exists(rds_path) &&
      !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_FB_FORCE"))) {
    message(sprintf("[FB-9.3.d-bootstrap] %s already present; skipping.",
                    entry$name))
    next
  }
  message(sprintf("[FB-9.3.d-bootstrap] fitting %s (family=%s K=%d p=%d seed=%d) ...",
                  entry$name, entry$family, entry$K, entry$p, entry$seed))
  t0 <- Sys.time()
  fit_result <- entry$fit_fn(entry)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  done in %.1f s; persisting RDS ...", elapsed))
  # Persist a draws-matrix snapshot mirroring the EB-side 8.6.D
  # convention (do NOT persist the full cmdstanr fit object because
  # it carries large temp file paths).
  draws_matrix <- posterior::as_draws_matrix(fit_result$fit$draws())
  snapshot <- list(
    golden_id        = entry$name,
    family           = entry$family,
    family_stan_id   = fit_result$family_stan_id,
    K                = entry$K,
    p                = entry$p,
    basis_type       = entry$basis_type,
    sub_phase        = "9.3.d",
    status           = "fitable",
    bootstrap_date   = bootstrap_date,
    cmdstan_version  = fit_result$cmdstan_version,
    R_version        = R_version_str,
    seed             = entry$seed,
    n_chains         = fit_result$chains,
    n_iter_warmup    = fit_result$iter_warmup,
    n_iter_sampling  = fit_result$iter_sampling,
    n_obs            = fit_result$n,
    draws_matrix     = draws_matrix,
    stan_data        = fit_result$stan_data,
    sampling_seconds = elapsed
  )
  saveRDS(snapshot, rds_path)
  gdpar:::.gdpar_golden_manifest_add(
    list(
      family          = entry$family,
      K               = entry$K,
      basis_type      = entry$basis_type,
      p               = entry$p,
      sub_phase       = "9.3.d",
      bootstrap_date  = bootstrap_date,
      cmdstan_version = cmdstan_version_str,
      R_version       = R_version_str,
      DHARMa_version  = NA_character_,
      fit_code_hash   = gdpar:::.gdpar_golden_fit_code_hash(entry$fit_fn),
      seed            = entry$seed,
      n_chains        = fit_result$chains,
      n_iter_warmup   = fit_result$iter_warmup,
      n_iter_sampling = fit_result$iter_sampling,
      rds_path        = rds_path_manifest,
      n_obs           = fit_result$n
    ),
    manifest_path
  )
  manifest_rows_written <- manifest_rows_written + 1L
}

message(sprintf(
  "[FB-9.3.d-bootstrap] DONE. Wrote %d FB KxP goldens + manifest rows under sub_phase = '9.3.d' (deuda D69 SUBSTANTIVELY CLOSED Sesion B9.6).",
  manifest_rows_written
))

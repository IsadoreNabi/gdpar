# Block RG, RG.2: R-native geometric sampling engine (Euclidean level).
#
# Default tests are pure R (no cmdstan): target/metric construction, the
# deterministic numerical layer (exact Hamiltonian flow on a Gaussian, energy
# conservation, reversibility, second-order accuracy) and the statistical layer
# (HMC recovers the moments of the easy suite targets). The cmdstan-backed
# layer (gradient cross-check closure-vs-$grad_log_prob, a cmdstan-backed
# target, the distributional match against NUTS) is gated by
# GDPAR_RUN_GEOMETRY_PILOTS, the same gate as the RG.1 pilots.

# --- helpers ---------------------------------------------------------------

# Exact Hamiltonian flow of the standard normal with the identity metric:
# dtheta/dt = p, dp/dt = -theta -> a rigid rotation by angle t in (theta, p).
exact_flow_std_normal <- function(theta0, p0, t) {
  list(theta = theta0 * cos(t) + p0 * sin(t),
       p = -theta0 * sin(t) + p0 * cos(t))
}

std_normal_target <- function(d) {
  gdpar_geom_target(log_prob = function(theta) -0.5 * sum(theta^2),
                    grad_log_prob = function(theta) -theta, dim = d)
}

# Neal's funnel (Block RG, RG.3): the canonical variable-curvature target whose
# remedy is a position-dependent (Riemannian) metric. The analytic expected
# Fisher information is block-diagonal -- 0.5(d-1)+1/s^2 for the log-scale v and
# exp(-v) for each of the d-1 nuisance coordinates -- and is the natural metric
# that resolves the neck. The analytic Hessian (of the log-density) drives the
# SoftAbs metric and the cmdstan $hessian cross-check.
funnel_target <- function(d = 10L, s = 3, with_hessian = FALSE) {
  lp <- function(th) {
    v <- th[1]; x <- th[-1]
    -0.5 * (v / s)^2 - 0.5 * (d - 1) * v - 0.5 * sum(x^2) * exp(-v)
  }
  gr <- function(th) {
    v <- th[1]; x <- th[-1]
    c(-v / s^2 - 0.5 * (d - 1) + 0.5 * sum(x^2) * exp(-v), -x * exp(-v))
  }
  he <- if (with_hessian) function(th) funnel_hessian(th, d, s) else NULL
  gdpar_geom_target(log_prob = lp, grad_log_prob = gr, dim = d, hessian = he)
}
funnel_fisher <- function(d = 10L, s = 3) {
  function(th) diag(c(0.5 * (d - 1) + 1 / s^2, rep(exp(-th[1]), d - 1)))
}
funnel_dfisher <- function(d = 10L, s = 3) {
  function(th) {
    L <- vector("list", d)
    L[[1]] <- diag(c(0, rep(-exp(-th[1]), d - 1)))
    for (k in 2:d) L[[k]] <- matrix(0, d, d)
    L
  }
}
funnel_hessian <- function(th, d = 10L, s = 3) {
  v <- th[1]; x <- th[-1]; idx <- 2:d
  H <- matrix(0, d, d)
  H[1, 1] <- -1 / s^2 - 0.5 * exp(-v) * sum(x^2)
  H[1, idx] <- x * exp(-v); H[idx, 1] <- x * exp(-v)
  H[cbind(idx, idx)] <- -exp(-v)
  H
}

# --- target construction ---------------------------------------------------

test_that("gdpar_geom_target builds from a closure and a suite instance", {
  tgt <- std_normal_target(2)
  expect_s3_class(tgt, "gdpar_geom_target")
  expect_identical(tgt$backend, "closure")
  expect_identical(tgt$dim, 2L)
  expect_equal(tgt$grad_log_prob(c(1, -2)), c(-1, 2))

  suite <- gdpar_geometry_suite()
  inst <- suite$G1_anisotropic$make(1, 100)
  tgt2 <- gdpar_geom_target(inst)
  expect_identical(tgt2$dim, inst$dim)
  expect_equal(tgt2$log_prob(rep(0, inst$dim)), inst$log_prob(rep(0, inst$dim)))
})

test_that("gdpar_geom_target validates its inputs", {
  expect_error(gdpar_geom_target(log_prob = function(x) 1, dim = 2),
               class = "gdpar_input_error")
  expect_error(
    gdpar_geom_target(log_prob = function(x) 1, grad_log_prob = function(x) x),
    class = "gdpar_input_error")
  expect_error(gdpar_geom_target(object = list(foo = 1)),
               class = "gdpar_input_error")
})

# --- Euclidean metric ------------------------------------------------------

test_that("gdpar_geom_metric_euclidean builds identity, diagonal and dense", {
  m <- gdpar_geom_metric_euclidean(dim = 3)
  expect_false(m$position_dependent)
  expect_equal(m$mass(rep(0, 3)), diag(3))
  expect_equal(m$inv_mass(rep(0, 3)), diag(3))
  expect_equal(m$logdet(rep(0, 3)), 0)

  md <- gdpar_geom_metric_euclidean(M = c(1, 4, 9))
  expect_equal(md$inv_mass(0) %*% md$mass(0), diag(3))
  expect_equal(md$logdet(0), sum(log(c(1, 4, 9))))
  # lower factor reproduces the mass: L L^T = M.
  L <- md$chol_mass(0)
  expect_equal(L %*% t(L), diag(c(1, 4, 9)))

  Sig <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  mD <- gdpar_geom_metric_euclidean(M = Sig)
  expect_equal(mD$inv_mass(0) %*% Sig, diag(2))
})

test_that("gdpar_geom_metric_euclidean rejects non-SPD and bad input", {
  expect_error(gdpar_geom_metric_euclidean(M = c(1, -1)),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_metric_euclidean(M = matrix(c(1, 2, 2, 1), 2)),
               class = "gdpar_input_error")  # indefinite
  expect_error(gdpar_geom_metric_euclidean(), class = "gdpar_input_error")
})

# --- numerical layer: deterministic, seed-free -----------------------------

test_that("leapfrog matches the exact Gaussian flow to second order", {
  tgt <- std_normal_target(3)
  metric <- gdpar_geom_metric_euclidean(dim = 3)
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(metric)
  theta0 <- c(0.7, -1.1, 0.3); p0 <- c(-0.2, 0.5, 1.0)
  Tend <- 1.2
  err <- function(eps) {
    L <- round(Tend / eps)
    end <- gdpar:::.gdpar_geom_leapfrog_traj(theta0, p0, tgt, metric, kin,
                                             eps, L)
    ex <- exact_flow_std_normal(theta0, p0, eps * L)
    sqrt(sum((end$theta - ex$theta)^2) + sum((end$p - ex$p)^2))
  }
  e1 <- err(0.04); e2 <- err(0.02)
  # halving the step quarters the error (order 2); allow a tolerance band.
  expect_gt(e1 / e2, 3.2)
  expect_lt(e1 / e2, 4.8)
  expect_lt(e2, 1e-2)
})

test_that("leapfrog conserves energy without drift (bounded, O(eps^2))", {
  tgt <- std_normal_target(2)
  metric <- gdpar_geom_metric_euclidean(dim = 2)
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(metric)
  H <- function(th, p) gdpar:::.gdpar_geom_hamiltonian(tgt, kin, th, p)
  max_dev <- function(eps, L) {
    th <- c(1.0, -0.5); p <- c(0.3, 0.8)
    H0 <- H(th, p); dev <- 0
    for (i in seq_len(L)) {
      st <- gdpar:::.gdpar_geom_leapfrog_step(th, p, tgt, metric, kin, eps)
      th <- st$theta; p <- st$p
      dev <- max(dev, abs(H(th, p) - H0))
    }
    dev
  }
  d1 <- max_dev(0.1, 400); d2 <- max_dev(0.05, 800)
  expect_lt(d2, d1)               # smaller step, smaller deviation
  expect_gt(d1 / d2, 2.5)         # ~4x for O(eps^2)
  expect_lt(d1, 0.5)              # bounded, no secular drift
})

test_that("leapfrog is exactly reversible under momentum flip", {
  suite <- gdpar_geometry_suite()
  inst <- suite$G1_anisotropic$make(1, 50)
  tgt <- gdpar_geom_target(inst)
  Sig <- diag(exp(seq(0, log(50), length.out = inst$dim)))
  metric <- gdpar_geom_metric_euclidean(M = Sig)
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(metric)
  set.seed(11)
  th0 <- stats::rnorm(inst$dim); p0 <- stats::rnorm(inst$dim)
  fwd <- gdpar:::.gdpar_geom_leapfrog_traj(th0, p0, tgt, metric, kin, 0.05, 30)
  back <- gdpar:::.gdpar_geom_leapfrog_traj(fwd$theta, -fwd$p, tgt, metric, kin,
                                            0.05, 30)
  expect_lt(max(abs(back$theta - th0)), 1e-9)
  expect_lt(max(abs(back$p + p0)), 1e-9)
})

# --- RG.3 Riemannian metric: construction ----------------------------------

test_that("SoftAbs Riemannian metric reduces to the identity on a normal", {
  tgt <- std_normal_target(3)
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "softabs", alpha = 1e6)
  expect_true(m$position_dependent)
  M <- m$mass(c(0.3, -0.4, 0.8))
  expect_lt(max(abs(M - diag(3))), 1e-4)              # softabs(1) ~ 1
  expect_equal(m$inv_mass(c(0.3, -0.4, 0.8)) %*% M, diag(3), tolerance = 1e-4)
  expect_lt(abs(m$logdet(c(0.1, 0.2, 0.3))), 1e-4)
})

test_that("Fisher Riemannian metric matches the analytic funnel information", {
  d <- 6L; s <- 2
  tgt <- funnel_target(d, s)
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "fisher",
                                    fisher = funnel_fisher(d, s),
                                    dfisher = funnel_dfisher(d, s))
  th <- c(-0.7, stats::rnorm(d - 1))
  expect_equal(m$mass(th), funnel_fisher(d, s)(th))
  # logdet of the diagonal information matrix.
  expect_equal(m$logdet(th),
               log(0.5 * (d - 1) + 1 / s^2) + (d - 1) * (-th[1]))
})

# --- RG.3 algebraic layer: metric derivatives vs finite differences --------

test_that("SoftAbs metric derivative matches finite differences", {
  d <- 5L; s <- 2
  tgt <- funnel_target(d, s, with_hessian = TRUE)   # exact Hessian -> clean dM
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "softabs", alpha = 50)
  th <- c(-0.5, 0.4, -0.2, 0.7, 0.1)[seq_len(d)]
  dM <- m$dmass(th)
  dM_fd <- gdpar:::.gdpar_geom_fd_dmass(m$mass, th, 1e-5)
  expect_lt(max(abs(unlist(dM) - unlist(dM_fd))), 1e-4)
})

test_that("Fisher metric derivative matches finite differences", {
  d <- 7L; s <- 3
  tgt <- funnel_target(d, s)
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "fisher",
                                    fisher = funnel_fisher(d, s),
                                    dfisher = funnel_dfisher(d, s))
  th <- c(-1.1, stats::rnorm(d - 1, 0, 0.5))
  dM <- m$dmass(th)
  dM_fd <- gdpar:::.gdpar_geom_fd_dmass(m$mass, th, 1e-5)
  expect_lt(max(abs(unlist(dM) - unlist(dM_fd))), 1e-7)
})

test_that("position-dependent kinetic gradient matches finite differences", {
  d <- 6L; s <- 2
  tgt <- funnel_target(d, s)
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "fisher",
                                    fisher = funnel_fisher(d, s),
                                    dfisher = funnel_dfisher(d, s))
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(m)
  set.seed(3); th <- c(-0.6, stats::rnorm(d - 1)); p <- stats::rnorm(d)
  g <- kin$grad_theta(th, p)
  h <- 1e-6
  g_fd <- vapply(seq_len(d), function(i) {
    e <- numeric(d); e[i] <- h
    (kin$value(th + e, p) - kin$value(th - e, p)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(g - g_fd)), 1e-5)
})

# --- RG.3 numerical layer: implicit generalised leapfrog -------------------

test_that("implicit generalised leapfrog is exactly reversible", {
  d <- 8L; s <- 3
  tgt <- funnel_target(d, s)
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "fisher",
                                    fisher = funnel_fisher(d, s),
                                    dfisher = funnel_dfisher(d, s))
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(m)
  set.seed(5); th0 <- c(-1.0, stats::rnorm(d - 1, 0, exp(-0.5)))
  p0 <- kin$draw_momentum(th0)
  fwd <- gdpar:::.gdpar_geom_leapfrog_step(th0, p0, tgt, m, kin, 0.2,
                                           fp_tol = 1e-12, fp_max = 200L)
  back <- gdpar:::.gdpar_geom_leapfrog_step(fwd$theta, -fwd$p, tgt, m, kin, 0.2,
                                            fp_tol = 1e-12, fp_max = 200L)
  expect_true(fwd$converged && back$converged)
  expect_lt(max(abs(back$theta - th0)), 1e-9)
  expect_lt(max(abs(back$p + p0)), 1e-9)
})

test_that("implicit generalised leapfrog conserves energy without drift", {
  d <- 8L; s <- 3
  tgt <- funnel_target(d, s)
  m <- gdpar_geom_metric_riemannian(tgt, curvature = "fisher",
                                    fisher = funnel_fisher(d, s),
                                    dfisher = funnel_dfisher(d, s))
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(m)
  H <- function(th, p) gdpar:::.gdpar_geom_hamiltonian(tgt, kin, th, p)
  set.seed(8); th0 <- c(-1.2, stats::rnorm(d - 1, 0, exp(-0.6)))
  p0 <- kin$draw_momentum(th0)
  max_dev <- function(eps, L) {
    th <- th0; p <- p0; H0 <- H(th, p); dev <- 0
    for (i in seq_len(L)) {
      st <- gdpar:::.gdpar_geom_leapfrog_step(th, p, tgt, m, kin, eps,
                                              fp_tol = 1e-12, fp_max = 200L)
      th <- st$theta; p <- st$p
      dev <- max(dev, abs(H(th, p) - H0))
    }
    dev
  }
  d1 <- max_dev(0.1, 60); d2 <- max_dev(0.05, 120)
  expect_lt(d1, 0.1)        # bounded, no secular drift
  expect_lt(d2, d1)         # smaller step, smaller deviation
})

# --- statistical layer: HMC recovers the easy targets ----------------------

test_that("HMC recovers the standard normal moments", {
  tgt <- std_normal_target(3)
  fit <- gdpar_geom_hmc(tgt, epsilon = 0.25, L = 18, n_iter = 4000,
                        n_warmup = 1000, seed = 123)
  expect_s3_class(fit, "gdpar_geom_hmc")
  expect_gt(fit$accept_rate, 0.7)
  expect_lt(max(abs(colMeans(fit$draws))), 0.12)
  expect_lt(max(abs(apply(fit$draws, 2, stats::var) - 1)), 0.2)
})

test_that("HMC recovers an anisotropic target; a matched dense metric helps", {
  suite <- gdpar_geometry_suite()
  inst <- suite$G1_anisotropic$make(1, 100)
  tgt <- gdpar_geom_target(inst)
  vars <- exp(seq(0, log(100), length.out = inst$dim))
  # Dense metric matched to the (known) covariance: efficient, high acceptance.
  matched <- gdpar_geom_metric_euclidean(M = diag(1 / vars))
  fit <- gdpar_geom_hmc(tgt, metric = matched, epsilon = 0.6, L = 20,
                        n_iter = 4000, n_warmup = 1000, seed = 7)
  expect_gt(fit$accept_rate, 0.6)
  # recover the marginal variances (the anisotropy ladder).
  rec <- apply(fit$draws, 2, stats::var)
  expect_lt(max(abs(log(rec / vars))), 0.6)
})

test_that("gdpar_geom_hmc validates its inputs", {
  tgt <- std_normal_target(2)
  expect_error(gdpar_geom_hmc(tgt, epsilon = -1), class = "gdpar_input_error")
  expect_error(gdpar_geom_hmc(tgt, n_warmup = -1), class = "gdpar_input_error")
  expect_error(gdpar_geom_hmc(tgt, init = c(0, 0, 0)),
               class = "gdpar_input_error")
})

# --- cmdstan-backed layer (gated) ------------------------------------------

test_that("[gated] closure gradients match cmdstan $grad_log_prob", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable cmdstan cross-checks.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  suite <- gdpar_geometry_suite()
  set.seed(404)
  for (nm in c("G0_isotropic", "G1_anisotropic", "G2_funnel")) {
    inst <- suite[[nm]]$make(suite[[nm]]$default_n,
                             suite[[nm]]$default_difficulty)
    f <- cmdstanr::write_stan_file(inst$stan_code)
    mod <- cmdstanr::cmdstan_model(f, compile_model_methods = TRUE,
                                   force_recompile = FALSE)
    fit <- mod$sample(data = inst$stan_data, chains = 1, iter_warmup = 1,
                      iter_sampling = 1, refresh = 0, show_messages = FALSE,
                      show_exceptions = FALSE)
    for (k in 1:3) {
      theta <- stats::rnorm(inst$dim, sd = 0.4)
      g_closure <- inst$grad_log_prob(theta)
      g_stan <- as.numeric(fit$grad_log_prob(unconstrained_variables = theta))
      expect_lt(max(abs(g_closure - g_stan)), 1e-6)
    }
  }
})

test_that("[gated] a cmdstan-backed target drives the same R leapfrog", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable cmdstan cross-checks.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  suite <- gdpar_geometry_suite()
  inst <- suite$G0_isotropic$make(1, 4)
  f <- cmdstanr::write_stan_file(inst$stan_code)
  mod <- cmdstanr::cmdstan_model(f, compile_model_methods = TRUE,
                                 force_recompile = TRUE)
  fit <- mod$sample(data = inst$stan_data, chains = 1, iter_warmup = 1,
                    iter_sampling = 1, refresh = 0, show_messages = FALSE,
                    show_exceptions = FALSE)
  tgt_stan <- gdpar_geom_target(fit, dim = inst$dim)
  tgt_closure <- gdpar_geom_target(inst)
  metric <- gdpar_geom_metric_euclidean(dim = inst$dim)
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(metric)
  th0 <- c(0.3, -0.6, 0.1, 0.9); p0 <- c(0.2, 0.4, -0.1, 0.5)
  a <- gdpar:::.gdpar_geom_leapfrog_traj(th0, p0, tgt_stan, metric, kin,
                                         0.1, 10)
  b <- gdpar:::.gdpar_geom_leapfrog_traj(th0, p0, tgt_closure, metric, kin,
                                         0.1, 10)
  expect_lt(max(abs(a$theta - b$theta)), 1e-6)
  expect_lt(max(abs(a$p - b$p)), 1e-6)
})

test_that("[gated] R HMC agrees distributionally with cmdstan NUTS on G0", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable cmdstan cross-checks.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  suite <- gdpar_geometry_suite()
  inst <- suite$G0_isotropic$make(1, 3)
  tgt <- gdpar_geom_target(inst)
  rfit <- gdpar_geom_hmc(tgt, epsilon = 0.3, L = 15, n_iter = 4000,
                         n_warmup = 1000, seed = 99)
  f <- cmdstanr::write_stan_file(inst$stan_code)
  mod <- cmdstanr::cmdstan_model(f)
  nuts <- mod$sample(data = inst$stan_data, chains = 2, iter_warmup = 500,
                     iter_sampling = 2000, refresh = 0, seed = 99,
                     show_messages = FALSE, show_exceptions = FALSE)
  nd <- posterior::as_draws_matrix(nuts$draws(variables = "theta"))
  expect_lt(max(abs(colMeans(rfit$draws) - colMeans(nd))), 0.15)
  expect_lt(max(abs(apply(rfit$draws, 2, stats::sd) -
                    apply(nd, 2, stats::sd))), 0.15)
})

# --- RG.3 statistical layer (gated, heavier): the funnel remedy ------------

test_that("[gated] RMHMC explores Neal's funnel where Euclidean HMC stalls", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the heavier funnel run.")
  }
  d <- 10L; s <- 3                       # true sd(v) = 3.
  tgt <- funnel_target(d, s)
  mF <- gdpar_geom_metric_riemannian(tgt, curvature = "fisher",
                                     fisher = funnel_fisher(d, s),
                                     dfisher = funnel_dfisher(d, s))
  eucl <- gdpar_geom_hmc(tgt, epsilon = 0.1, L = 20, n_iter = 800,
                         n_warmup = 400, seed = 42)
  rmhmc <- gdpar_geom_hmc(tgt, metric = mF, epsilon = 0.3, L = 12,
                          n_iter = 800, n_warmup = 400, seed = 42)
  # The energy diagnostic exposes the Euclidean metric as ill-matched (E-BFMI
  # far below the 0.3 rule of thumb) while the Riemannian metric resolves it.
  expect_lt(eucl$ebfmi, 0.2)
  expect_gt(rmhmc$ebfmi, 0.5)
  expect_gt(rmhmc$ebfmi, 3 * eucl$ebfmi)
  # The Riemannian sampler reaches deep into the neck and recovers Var(v) = s^2.
  expect_lt(min(rmhmc$draws[, 1]), min(eucl$draws[, 1]) - 1)
  expect_lt(abs(stats::sd(rmhmc$draws[, 1]) - s),
            abs(stats::sd(eucl$draws[, 1]) - s))
  expect_lt(abs(stats::sd(rmhmc$draws[, 1]) - s), 0.4)
})

test_that("[gated] cmdstan $hessian matches the analytic Hessian of the suite", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable cmdstan cross-checks.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  suite <- gdpar_geometry_suite()
  analytic <- list(
    G0_isotropic = function(inst) function(th) -diag(inst$dim),
    G1_anisotropic = function(inst) {
      iv <- 1 / exp(seq(0, log(100), length.out = inst$dim))
      function(th) -diag(iv)
    },
    G2_funnel = function(inst) function(th) funnel_hessian(th, inst$dim, 3))
  set.seed(2024)
  for (nm in names(analytic)) {
    inst <- suite[[nm]]$make(suite[[nm]]$default_n, suite[[nm]]$default_difficulty)
    f <- cmdstanr::write_stan_file(inst$stan_code)
    # force_recompile guarantees the model methods ($hessian) are built in even
    # if another test left a methods-less executable in the cache.
    mod <- cmdstanr::cmdstan_model(f, compile_model_methods = TRUE,
                                   force_recompile = TRUE)
    fit <- mod$sample(data = inst$stan_data, chains = 1, iter_warmup = 1,
                      iter_sampling = 1, refresh = 0, show_messages = FALSE,
                      show_exceptions = FALSE)
    hfun <- analytic[[nm]](inst)
    for (k in 1:3) {
      theta <- stats::rnorm(inst$dim, sd = 0.3)
      h_stan <- fit$hessian(unconstrained_variables = theta)$hessian
      expect_lt(max(abs(h_stan - hfun(theta))), 1e-4)
    }
    # SoftAbs metric built from the Stan Hessian equals the analytic-Hessian one.
    tgt_stan <- gdpar_geom_target(fit, dim = inst$dim)
    tgt_an <- gdpar_geom_target(log_prob = inst$log_prob,
                                grad_log_prob = inst$grad_log_prob,
                                dim = inst$dim, hessian = hfun)
    mS_stan <- gdpar_geom_metric_riemannian(tgt_stan, curvature = "softabs")
    mS_an <- gdpar_geom_metric_riemannian(tgt_an, curvature = "softabs")
    th <- stats::rnorm(inst$dim, sd = 0.3)
    expect_lt(max(abs(mS_stan$mass(th) - mS_an$mass(th))), 1e-4)
  }
})

# --- RG.3 Capa 2: learned Gaussian-process expected-Fisher metric ----------

# Algebraic layer A: the log-Cholesky parametrisation round-trips and both its
# forward (dM from dpsi) and inverse (dpsi from dM) differentials are exact.
test_that("log-Cholesky parametrisation round-trips and differentiates", {
  set.seed(1)
  d <- 4L
  idx <- gdpar:::.gdpar_geom_tri_index(d)
  A <- matrix(stats::rnorm(d * d), d); M <- crossprod(A) + diag(d)
  lc <- gdpar:::.gdpar_geom_logchol(M, idx)
  expect_lt(max(abs(gdpar:::.gdpar_geom_L_from_psi(lc$psi, idx, d) - lc$L)), 1e-12)
  expect_lt(max(abs(tcrossprod(lc$L) - M)), 1e-10)
  # dM_from_dpsi vs finite differences of M(psi) along a random direction.
  dpsi <- stats::rnorm(nrow(idx)); hh <- 1e-6
  Lp <- gdpar:::.gdpar_geom_L_from_psi(lc$psi + hh * dpsi, idx, d)
  Lm <- gdpar:::.gdpar_geom_L_from_psi(lc$psi - hh * dpsi, idx, d)
  dM_fd <- (tcrossprod(Lp) - tcrossprod(Lm)) / (2 * hh)
  dM <- gdpar:::.gdpar_geom_dM_from_dpsi(lc$L, dpsi, idx, d)
  expect_lt(max(abs(dM - dM_fd)), 1e-6)
  # dpsi_from_dM is the inverse of dM_from_dpsi.
  expect_lt(max(abs(gdpar:::.gdpar_geom_dpsi_from_dM(lc$L, dM, idx, d) - dpsi)),
            1e-8)
})

# Algebraic layer B: with a small nugget the surrogate interpolates the supplied
# Fisher at the reservoir sites, and the mass is SPD everywhere by construction.
test_that("the learned metric interpolates the Fisher at the reservoir sites", {
  set.seed(2); d <- 3L
  tgt <- std_normal_target(d)
  fish <- function(theta) diag(c(1 + 0.5 * theta[1]^2, 2, 0.5))
  sites <- matrix(stats::rnorm(40 * d, sd = 0.8), ncol = d)
  m <- gdpar_geom_metric_gp_fisher(tgt, fisher = fish, sites = sites,
                                   nugget = 1e-8)
  expect_true(m$position_dependent)
  expect_identical(m$metric_kind, "gp_fisher")
  expect_identical(m$n_sites, 40L)
  s <- sites[7, ]
  expect_lt(max(abs(m$mass(s) - fish(s))), 1e-3)            # interpolation.
  ev <- eigen(m$mass(c(0.5, -0.3, 0.2)), symmetric = TRUE,
              only.values = TRUE)$values
  expect_gt(min(ev), 0)                                    # SPD by construction.
})

# Algebraic layer C: the closed-form metric derivative (analytic kernel
# derivative + Daleckii-Krein SoftAbs derivative through the log-Cholesky map)
# matches finite differences of the mass.
test_that("the learned metric derivative matches finite differences", {
  set.seed(3); d <- 4L; s <- 2
  tgt <- funnel_target(d, s, with_hessian = TRUE)
  sites <- cbind(stats::runif(30, -3, 1),
                 matrix(stats::rnorm(30 * (d - 1), 0, 0.5), ncol = d - 1))
  m <- gdpar_geom_metric_gp_fisher(tgt, fisher = funnel_fisher(d, s),
                                   sites = sites, alpha = 50, nugget = 1e-6)
  th <- c(-0.6, 0.3, -0.2, 0.5)
  dM <- m$dmass(th)
  dM_fd <- gdpar:::.gdpar_geom_fd_dmass(m$mass, th, 1e-5)
  expect_lt(max(abs(unlist(dM) - unlist(dM_fd))), 1e-3)
})

# Methodological layer: graceful continuous degradation to SoftAbs far from the
# reservoir, and the predictive standard deviation as the novelty detector.
test_that("the learned metric degrades to SoftAbs far away and flags novelty", {
  set.seed(4); d <- 3L
  tgt <- std_normal_target(d)
  fish <- function(theta) diag(rep(3, d))   # clearly distinct from SoftAbs (= I).
  sites <- matrix(stats::rnorm(30 * d, sd = 0.5), ncol = d)
  m <- gdpar_geom_metric_gp_fisher(tgt, fisher = fish, sites = sites,
                                   nugget = 1e-8)
  mS <- gdpar_geom_metric_riemannian(tgt, curvature = "softabs")
  near <- sites[1, ]; far <- rep(50, d)
  expect_lt(max(abs(m$mass(near) - fish(near))), 1e-2)     # follows the Fisher.
  expect_lt(max(abs(m$mass(far) - mS$mass(far))), 1e-6)    # returns to SoftAbs.
  expect_gt(m$novelty(far), 0.99)
  expect_lt(m$novelty(near), 0.05)
})

# Numerical layer: the implicit generalised leapfrog with the learned metric is
# exactly reversible, and conserves the Hamiltonian without secular drift.
test_that("implicit leapfrog with the learned metric is reversible and conserves H", {
  set.seed(5); d <- 3L; s <- 2
  # Reversibility on the funnel (the hard, position-dependent case): exact up to
  # the fixed-point tolerance wherever the implicit solve converges. The GP
  # metric (SoftAbs mean + learned residual) is stiffer than the analytic
  # diagonal Fisher, so the fixed point converges at a smaller step than Capa 1.
  tgt <- funnel_target(d, s, with_hessian = TRUE)
  sites <- cbind(stats::runif(25, -3, 1),
                 matrix(stats::rnorm(25 * (d - 1), 0, 0.5), ncol = d - 1))
  m <- gdpar_geom_metric_gp_fisher(tgt, fisher = funnel_fisher(d, s),
                                   sites = sites, alpha = 50, nugget = 1e-6)
  kin <- gdpar:::.gdpar_geom_kinetic_gaussian(m)
  th0 <- c(-0.5, 0.2, -0.1); p0 <- kin$draw_momentum(th0)
  fwd <- gdpar:::.gdpar_geom_leapfrog_step(th0, p0, tgt, m, kin, 0.05,
                                           fp_tol = 1e-12, fp_max = 200L)
  back <- gdpar:::.gdpar_geom_leapfrog_step(fwd$theta, -fwd$p, tgt, m, kin, 0.05,
                                            fp_tol = 1e-12, fp_max = 200L)
  expect_true(fwd$converged && back$converged)
  expect_lt(max(abs(back$theta - th0)), 1e-8)
  expect_lt(max(abs(back$p + p0)), 1e-8)
  # Energy conservation without secular drift on a benign, well-covered smooth
  # learned metric. A non-converged step is rejected by the sampler, so we
  # accumulate the deviation only over the converged trajectory prefix.
  tgtN <- std_normal_target(d)
  fishN <- function(theta) diag(c(1 + 0.5 * theta[1]^2, 2, 0.5))
  mN <- gdpar_geom_metric_gp_fisher(
    tgtN, fisher = fishN,
    sites = matrix(stats::rnorm(40 * d, sd = 0.8), ncol = d), nugget = 1e-8)
  kinN <- gdpar:::.gdpar_geom_kinetic_gaussian(mN)
  HN <- function(th, p) gdpar:::.gdpar_geom_hamiltonian(tgtN, kinN, th, p)
  th <- c(0.2, -0.3, 0.1); p <- kinN$draw_momentum(th)
  H0 <- HN(th, p); dev <- 0; nconv <- 0L
  for (i in seq_len(30)) {
    stp <- gdpar:::.gdpar_geom_leapfrog_step(th, p, tgtN, mN, kinN, 0.1,
                                             fp_tol = 1e-12, fp_max = 200L)
    if (!isTRUE(stp$converged)) break
    th <- stp$theta; p <- stp$p
    dev <- max(dev, abs(HN(th, p) - H0)); nconv <- nconv + 1L
  }
  expect_gt(nconv, 20L)        # the integrator runs the full trajectory,
  expect_lt(dev, 0.05)         # conserving energy without drift.
})

test_that("gdpar_geom_metric_gp_fisher validates its inputs", {
  tgt <- std_normal_target(2)
  expect_error(gdpar_geom_metric_gp_fisher(tgt, fisher = 1,
                                           sites = matrix(0, 2, 2)),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_metric_gp_fisher(tgt, fisher = function(th) diag(2),
                                           sites = matrix(0, 2, 3)),
               class = "gdpar_input_error")   # wrong number of columns.
  expect_error(gdpar_geom_metric_gp_fisher(tgt, fisher = function(th) diag(2),
                                           sites = matrix(0, 1, 2)),
               class = "gdpar_input_error")   # too few sites.
})

test_that("gdpar_geom_reservoir collects positions from a warmup run", {
  tgt <- std_normal_target(3)
  sites <- gdpar_geom_reservoir(tgt, n_sites = 20, n_warmup = 80, seed = 7)
  expect_equal(dim(sites), c(20L, 3L))
  expect_true(all(is.finite(sites)))
})

# Statistical layer (gated, heavier): the learned metric reproduces the funnel
# exploration that the analytic Fisher achieved in Capa 1.
test_that("[gated] the learned GP metric explores Neal's funnel like the analytic Fisher", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the heavier learned-funnel run.")
  }
  d <- 10L; s <- 3                              # true sd(v) = 3.
  tgt <- funnel_target(d, s, with_hessian = TRUE)
  # Phase 1: a reservoir spanning the log-scale v (the funnel's curvature axis).
  set.seed(202)
  vgrid <- rep(seq(-6, 2, length.out = 20), each = 3)
  sites <- cbind(vgrid,
                 matrix(stats::rnorm(length(vgrid) * (d - 1), 0, 0.5),
                        ncol = d - 1))
  # A moderate SoftAbs softening keeps the mean (hence the metric and its
  # derivative) smooth enough for the implicit solve; the learned residual then
  # supplies the funnel's expected-Fisher curvature.
  m <- gdpar_geom_metric_gp_fisher(tgt, fisher = funnel_fisher(d, s),
                                   sites = sites, alpha = 50, nugget = 1e-6)
  eucl <- gdpar_geom_hmc(tgt, epsilon = 0.1, L = 20, n_iter = 600,
                         n_warmup = 300, seed = 42)
  # the GP metric is stiffer than the analytic diagonal Fisher, so the implicit
  # fixed point converges through the neck at a smaller step (longer trajectory).
  rmhmc <- gdpar_geom_hmc(tgt, metric = m, epsilon = 0.02, L = 40,
                          n_iter = 600, n_warmup = 300, seed = 42)
  # The decisive, robustly achievable signals (matching Capa 1's framing): the
  # energy diagnostic and the depth of neck exploration. The euclidean metric
  # stalls (E-BFMI far below 0.3); the learned metric lifts it towards 1 and
  # reaches deep into the neck where the euclidean cannot. Full Var(v) recovery
  # needs far longer runs at this small step and is not asserted here.
  expect_gt(rmhmc$accept_rate, 0.5)
  expect_lt(eucl$ebfmi, 0.2)
  expect_gt(rmhmc$ebfmi, 0.5)
  expect_gt(rmhmc$ebfmi, 2 * eucl$ebfmi)
  expect_lt(min(rmhmc$draws[, 1]), min(eucl$draws[, 1]) - 1)
})

# ===========================================================================
# RG.3 Capa 2 (completion): simulation-based expected Fisher + active learning.
#
# Algebraic layer (ungated, deterministic by seed): the score-outer-product
# estimator is positive semi-definite and recovers the closed-form expected
# Fisher -- a constant one (a normal location model) and a position-dependent
# one (a Poisson log-rate model). Statistical layer for the online active
# learning loop (ungated tiny smoke; gated heavier recovery + a real
# cmdstan-backed count model, the rehearsal for the eBird Tweedie of RG.7).
# ===========================================================================

# --- generative helpers ----------------------------------------------------

# Multivariate normal location model y ~ N(theta, Sigma0): expected Fisher is
# the constant precision solve(Sigma0).
gaussian_location_target <- function(Sigma0) {
  P0 <- solve(Sigma0); Lc <- t(chol(Sigma0)); d <- nrow(Sigma0)
  gdpar_geom_target(
    log_prob = function(theta) -0.5 * as.numeric(crossprod(theta, P0 %*% theta)),
    grad_log_prob = function(theta) -as.numeric(P0 %*% theta),
    hessian = function(theta) -P0, dim = d,
    simulate = function(theta) as.numeric(theta + Lc %*% stats::rnorm(d)),
    score = function(theta, y) as.numeric(P0 %*% (y - theta)))
}

# Poisson log-rate with n_obs observations: loglik score = sum(y) - n exp(theta),
# expected Fisher = n exp(theta) (grows with theta -- position-dependent).
poisson_rate_target <- function(n_obs = 25L) {
  gdpar_geom_target(
    log_prob = function(theta) n_obs * (3 * theta - exp(theta)),
    grad_log_prob = function(theta) n_obs * (3 - exp(theta)), dim = 1,
    simulate = function(theta) stats::rpois(n_obs, exp(theta)),
    score = function(theta, y) sum(y) - n_obs * exp(theta))
}

# --- algebraic layer: the simulation-based expected Fisher estimator --------

test_that("the score-outer-product estimator is PSD and recovers a constant Fisher", {
  Sigma0 <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
  tgt <- gaussian_location_target(Sigma0)
  fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 6000, seed = 1)
  Ihat <- fisher(c(0.3, -0.2))
  expect_equal(attr(Ihat, "n_sim"), 6000L)
  # symmetric, positive-definite, and close to the closed-form precision.
  expect_lt(max(abs(Ihat - t(Ihat))), 1e-12)
  expect_gt(min(eigen(Ihat, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_lt(max(abs(Ihat - solve(Sigma0))), 0.08)
  # deterministic function of theta (reproducible regardless of call order).
  expect_identical(fisher(c(0.3, -0.2)), fisher(c(0.3, -0.2)))
})

test_that("the estimator recovers a position-dependent Fisher (Poisson log-rate)", {
  n_obs <- 25L
  tgt <- poisson_rate_target(n_obs)
  fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 6000, seed = 7)
  ests <- vapply(c(-0.5, 0, 0.5), function(th) as.numeric(fisher(th)), numeric(1))
  truth <- n_obs * exp(c(-0.5, 0, 0.5))
  expect_lt(max(abs(ests - truth) / truth), 0.06)   # within 6% at n_sim = 6000.
  expect_true(all(diff(ests) > 0))                  # grows with theta.
})

test_that("a rank-deficient simulation (n_sim < dim) is floored to SPD", {
  Sigma0 <- diag(c(1, 2, 3))
  tgt <- gaussian_location_target(Sigma0)
  fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 1L, seed = 3, floor = 1e-8)
  Ihat <- fisher(c(0, 0, 0))
  expect_gte(min(eigen(Ihat, symmetric = TRUE, only.values = TRUE)$values), 1e-9)
})

test_that("gdpar_geom_fisher_simulator validates its inputs", {
  ng <- gdpar_geom_target(log_prob = function(t) -0.5 * sum(t^2),
                          grad_log_prob = function(t) -t, dim = 2)
  expect_error(gdpar_geom_fisher_simulator(ng),
               class = "gdpar_input_error")          # not generative.
  bad <- gdpar_geom_target(log_prob = function(t) -0.5 * sum(t^2),
                           grad_log_prob = function(t) -t, dim = 2,
                           simulate = function(theta) stats::rnorm(2),
                           score = function(theta, y) c(y, y))  # wrong length.
  f <- gdpar_geom_fisher_simulator(bad, n_sim = 2)
  expect_error(f(c(0, 0)), class = "gdpar_input_error")
})

test_that("the simulated Fisher feeds the GP surrogate (mass SPD, novelty far)", {
  Sigma0 <- matrix(c(1, 0.4, 0.4, 1.5), 2, 2)
  tgt <- gaussian_location_target(Sigma0)
  fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 3000, seed = 5)
  sites <- gdpar_geom_reservoir(tgt, n_sites = 18, n_warmup = 60, seed = 11)
  m <- gdpar_geom_metric_gp_fisher(tgt, fisher = fisher, sites = sites,
                                   nugget = 1e-3)
  M0 <- m$mass(c(0, 0))
  expect_gt(min(eigen(M0, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_lt(max(abs(M0 - solve(Sigma0))), 0.2)       # near data: close to truth.
  expect_lt(m$novelty(sites[1, ]), 0.2)              # at a site: low novelty.
  expect_gt(m$novelty(c(40, -40)), 0.9)              # far: high novelty.
})

# --- statistical layer: the online active-learning loop (ungated tiny) ------

test_that("gdpar_geom_rmhmc_adaptive runs, grows the reservoir, stays exact", {
  Sigma0 <- diag(c(1, 4))
  tgt <- gaussian_location_target(Sigma0)
  fit <- gdpar_geom_rmhmc_adaptive(
    tgt, n_sim = 300, n_sites_init = 10, max_rounds = 3, batch = 15,
    n_add = 6, epsilon = 0.25, L = 12, n_iter = 60, n_warmup = 30, seed = 42)
  expect_s3_class(fit, "gdpar_geom_rmhmc_adaptive")
  expect_true(all(is.finite(fit$draws)))
  expect_identical(fit$n_sites_trace[1], 10L)        # starts at n_sites_init.
  expect_true(all(diff(fit$n_sites_trace) >= 0))     # reservoir never shrinks.
  expect_length(fit$novelty_trace, fit$n_rounds)
  expect_true(fit$metric$position_dependent)         # a learned Riemannian metric.
  expect_gte(nrow(fit$reservoir), 10L)
  # the metric stays SPD on the explored region (exactness preconditioner).
  M <- fit$metric$mass(colMeans(fit$draws))
  expect_gt(min(eigen(M, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("gdpar_geom_rmhmc_adaptive validates its inputs", {
  tgt <- gaussian_location_target(diag(2))
  expect_error(gdpar_geom_rmhmc_adaptive(tgt, n_sites_init = 1L),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_rmhmc_adaptive(tgt, fisher = 3),
               class = "gdpar_input_error")
})

# --- statistical layer (gated, heavier): recovery + a real cmdstan count model

test_that("[gated] the adaptive loop recovers the moments of a normal location model", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the adaptive recovery run.")
  }
  Sigma0 <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
  tgt <- gaussian_location_target(Sigma0)
  fit <- gdpar_geom_rmhmc_adaptive(
    tgt, n_sim = 1500, n_sites_init = 30, max_rounds = 4, batch = 60,
    n_add = 20, epsilon = 0.2, L = 25, n_iter = 3000, n_warmup = 1000,
    seed = 7)
  expect_lt(max(abs(colMeans(fit$draws))), 0.15)            # mean (truth 0).
  expect_lt(max(abs(stats::cov(fit$draws) - Sigma0)), 0.3)  # covariance.
  expect_gt(fit$ebfmi, 0.3)
})

test_that("[gated] full pipeline on a real cmdstan-backed Poisson count model", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the cmdstan demo.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  # A genuine count likelihood with data -- the structural rehearsal for the
  # eBird Tweedie of RG.7. The sampler is driven by the REAL cmdstan
  # $grad_log_prob; the simulation-based Fisher uses the model's closed-form
  # likelihood score (validated against the analytic X^T W X below).
  set.seed(123)
  n <- 60L; x <- stats::rnorm(n)
  b_true <- c(0.4, 0.8)
  y <- stats::rpois(n, exp(b_true[1] + b_true[2] * x))
  X <- cbind(1, x)
  stan_code <- paste(
    "data { int<lower=1> n; vector[n] x; array[n] int<lower=0> y; }",
    "parameters { real b0; real b1; }",
    "model { b0 ~ normal(0, 5); b1 ~ normal(0, 5);",
    "        y ~ poisson_log(b0 + b1 * x); }", sep = "\n")
  f <- cmdstanr::write_stan_file(stan_code)
  mod <- cmdstanr::cmdstan_model(f, compile_model_methods = TRUE,
                                 force_recompile = TRUE)
  fit1 <- mod$sample(data = list(n = n, x = x, y = y), chains = 1,
                     iter_warmup = 1, iter_sampling = 1, refresh = 0,
                     show_messages = FALSE, show_exceptions = FALSE)
  eta <- function(b) X %*% b
  tgt <- gdpar_geom_target(
    fit1, dim = 2,
    simulate = function(b) stats::rpois(n, exp(as.numeric(eta(b)))),
    score = function(b, ysim)
      as.numeric(crossprod(X, ysim - exp(as.numeric(eta(b))))))
  expect_true(!is.null(tgt$simulate) && !is.null(tgt$score))
  # the simulated Fisher matches the analytic GLM information X^T W X.
  b0v <- c(0.4, 0.8)
  fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 8000, seed = 9)
  W <- diag(as.numeric(exp(eta(b0v))))
  expect_lt(max(abs(fisher(b0v) - crossprod(X, W %*% X))) /
              max(abs(crossprod(X, W %*% X))), 0.1)
  # end-to-end adaptive RMHMC over the cmdstan density.
  afit <- gdpar_geom_rmhmc_adaptive(
    tgt, n_sim = 1500, n_sites_init = 20, max_rounds = 3, batch = 40,
    n_add = 12, epsilon = 0.15, L = 20, n_iter = 1500, n_warmup = 600,
    init = c(0, 0), seed = 11)
  expect_true(all(is.finite(afit$draws)))
  expect_true(is.finite(afit$ebfmi))
  expect_lt(max(abs(colMeans(afit$draws) - b_true)), 0.25)   # recovers the truth.
})

# ===========================================================================
# RG.4: sub-Riemannian geometry -- the remedy for a quasi-deterministic
# posterior (the eBird count / tweedie case). The accessible distribution is
# read from the near-null space of the expected Fisher; a continuous spectral
# filter (no hard cut) splits floor from walls; the integrator is a Strang
# splitting whose stiff wall block flows EXACTLY (a closed-form harmonic
# rotation, no step-size penalty) while the gentle residual is leapfrogged and
# the Metropolis correction with the exact density keeps the sampler exact.
#
# Algebraic layer (ungated, deterministic): the exact flow conserves the
# reference quadratic energy to machine precision and is reversible; the Strang
# trajectory is reversible; the integrator is exact on a pure quadratic at any
# step; the verticality filter is monotone and bounded. Statistical layer:
# ungated (a mild canyon recovered through gdpar_geom_hmc) plus gated heavier
# runs (G4 at large n where Euclidean collapses, and a real cmdstan-backed
# quasi-deterministic count -- the rehearsal of the eBird tweedie of RG.7).
# ===========================================================================

# --- algebraic layer: the exact harmonic flow and the Strang integrator -----

test_that("the exact harmonic flow conserves the reference quadratic and reverses", {
  set.seed(1)
  d <- 4
  U <- qr.Q(qr(matrix(stats::rnorm(d * d), d, d)))   # a random orthonormal basis.
  omega <- c(0, 0.7, 5, 30)                          # a floor mode and stiff walls.
  ref <- c(0.2, -0.1, 0.4, 0.0)
  theta <- c(1, -2, 0.5, 1.3); p <- c(-0.3, 0.8, 1.1, -0.6)
  A <- U %*% (omega^2 * t(U))
  HQ <- function(th, pp)
    0.5 * sum(pp^2) + 0.5 * as.numeric(crossprod(th - ref, A %*% (th - ref)))
  st <- .gdpar_geom_subriemann_flow(theta, p, ref, U, omega, 0.37)
  # the closed-form flow conserves the harmonic energy to machine precision...
  expect_lt(abs(HQ(st$theta, st$p) - HQ(theta, p)), 1e-11)
  # ...and is exactly reversible: flow(-t) undoes flow(t).
  back <- .gdpar_geom_subriemann_flow(st$theta, st$p, ref, U, omega, -0.37)
  expect_lt(max(abs(back$theta - theta), abs(back$p - p)), 1e-11)
})

test_that("the Strang sub-Riemannian trajectory is time-reversible", {
  tgt <- gdpar_geom_target(
    log_prob = function(th) -0.5 * (th[1]^2 + 80 * th[2]^2) - 0.05 * th[1]^4,
    grad_log_prob = function(th)
      -c(th[1] + 0.2 * th[1]^3, 80 * th[2]), dim = 2)
  metric <- gdpar_geom_metric_subriemannian(
    tgt, fisher = function(th) diag(c(1, 80)))
  th0 <- c(0.6, 0.1); p0 <- c(-0.4, 0.25)
  fwd <- metric$integrator(th0, p0, tgt, 0.3, 8L)
  rev <- metric$integrator(fwd$theta, -fwd$p, tgt, 0.3, 8L)
  expect_true(fwd$converged && rev$converged)
  expect_lt(max(abs(rev$theta - th0)), 1e-9)      # position returns...
  expect_lt(max(abs(-rev$p - p0)), 1e-9)          # ...and so does the momentum.
})

test_that("the integrator is exact on a pure quadratic (all directions walls)", {
  P <- matrix(c(3, 1, 1, 200), 2, 2)              # a stiff Gaussian.
  tgt <- gdpar_geom_target(
    log_prob = function(th) -0.5 * as.numeric(crossprod(th, P %*% th)),
    grad_log_prob = function(th) -as.numeric(P %*% th), dim = 2)
  # force every direction to be a wall: tau far below the spectrum.
  metric <- gdpar_geom_metric_subriemannian(
    tgt, fisher = function(th) P, tau = 1e-10)
  H <- function(th, pp) 0.5 * sum(pp^2) -
    (-0.5 * as.numeric(crossprod(th, P %*% th)))
  th0 <- c(1.2, 0.3); p0 <- c(0.5, -0.7)
  prop <- metric$integrator(th0, p0, tgt, 1.5, 20L)   # a large step on purpose.
  expect_lt(abs(H(prop$theta, prop$p) - H(th0, p0)), 1e-6)  # exact: no drift.
})

test_that("the verticality filter is monotone, bounded, and caps the residual", {
  tgt <- gdpar_geom_target(log_prob = function(t) -0.5 * sum(t^2),
                           grad_log_prob = function(t) -t, dim = 3)
  lam <- c(1, 12, 900)
  metric <- gdpar_geom_metric_subriemannian(
    tgt, fisher = function(th) diag(lam))
  ord <- order(metric$eigenvalues)
  w <- metric$verticality[ord]; lv <- metric$eigenvalues[ord]
  expect_true(all(w > 0 & w < 1))                 # strictly interior (no hard cut).
  expect_true(all(diff(w) > 0))                   # monotone increasing in lambda.
  expect_identical(metric$metric_kind, "sub_riemannian")
  expect_false(metric$position_dependent)
  expect_equal(metric$n_walls, 2L)                # the two stiff directions.
  # residual curvature (1 - w) * lambda is capped at the floor scale (tau = 1).
  expect_lt(max((1 - metric$verticality) * metric$eigenvalues), 1 + 1e-8)
  expect_equal(metric$frequencies, sqrt(metric$verticality * metric$eigenvalues))
})

test_that("gdpar_geom_metric_subriemannian validates its inputs", {
  tgt <- gdpar_geom_target(log_prob = function(t) -0.5 * sum(t^2),
                           grad_log_prob = function(t) -t, dim = 2)
  expect_error(gdpar_geom_metric_subriemannian(tgt, fisher = 3),
               class = "gdpar_input_error")
  expect_error(
    gdpar_geom_metric_subriemannian(tgt, fisher = function(th) diag(2),
                                    reference = c(0, 0, 0)),
    class = "gdpar_input_error")
  expect_error(
    gdpar_geom_metric_subriemannian(tgt, fisher = function(th) diag(3)),
    class = "gdpar_input_error")                  # wrong-size Fisher.
})

# --- statistical layer (ungated, light): a mild canyon through gdpar_geom_hmc

test_that("the sub-Riemannian sampler recovers a mild canyon and exposes the path", {
  tgt <- gdpar_geom_target(
    log_prob = function(th) -0.5 * (th[1]^2 + 50 * th[2]^2),
    grad_log_prob = function(th) -c(th[1], 50 * th[2]), dim = 2)
  metric <- gdpar_geom_metric_subriemannian(
    tgt, fisher = function(th) diag(c(1, 50)))
  expect_true(is.function(metric$integrator))
  fit <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.5, L = 12,
                        n_iter = 2500, n_warmup = 600, seed = 3)
  expect_identical(fit$metric_type, "euclidean_constant")  # identity kinetic.
  expect_true(all(is.finite(fit$draws)))
  expect_gt(fit$accept_rate, 0.6)
  expect_equal(fit$n_divergent, 0L)
  expect_lt(abs(stats::sd(fit$draws[, 1]) - 1), 0.25)             # floor sd ~ 1.
  expect_lt(abs(stats::sd(fit$draws[, 2]) - 1 / sqrt(50)), 0.05)  # wall sd ~ 0.14.
})

# --- statistical layer (gated): G4 at large n -- Euclidean collapses ---------

test_that("[gated] the sub-Riemannian sampler solves the G4 canyon where Euclidean stalls", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the G4 canyon run.")
  }
  suite <- gdpar_geometry_suite("G4_quasi_deterministic")
  inst <- suite$G4_quasi_deterministic$make(n = 1000, difficulty = 4)
  tgt <- gdpar_geom_target(inst, dim = inst$dim)
  n_eff <- 1000
  fisher <- function(theta) diag(c(1, rep(n_eff, inst$dim - 1)))
  metric <- gdpar_geom_metric_subriemannian(tgt, fisher = fisher)
  expect_equal(metric$n_walls, inst$dim - 1L)
  # the canyon-aware sampler takes an O(1) step the stiff walls would forbid.
  sub <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.5, L = 12,
                        n_iter = 3000, n_warmup = 800, seed = 1)
  expect_equal(sub$n_divergent, 0L)
  expect_gt(sub$accept_rate, 0.7)
  expect_gt(sub$ebfmi, 0.3)
  expect_lt(abs(stats::sd(sub$draws[, 1]) - 1), 0.2)                # floor sd ~ 1.
  expect_lt(abs(stats::sd(sub$draws[, 2]) - 1 / sqrt(n_eff)) /
              (1 / sqrt(n_eff)), 0.4)                               # wall sd ~ n^-1/2.
  # at the same step the Euclidean sampler is wrecked by the walls (omega eps
  # far above the stability limit): it diverges and cannot move along the floor.
  euc <- gdpar_geom_hmc(tgt, epsilon = 0.5, L = 12, n_iter = 3000,
                        n_warmup = 800, seed = 1)
  expect_gt(euc$n_divergent, sub$n_divergent + 500L)
  expect_lt(stats::sd(euc$draws[, 1]), 0.5)                        # stuck on the floor.
})

# --- statistical layer (gated, real cmdstan): the quasi-deterministic count --

test_that("[gated] sub-Riemannian samples a real cmdstan-backed near-deterministic count", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the cmdstan count demo.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  # A near-deterministic Poisson count: a wide predictor makes exp(eta) huge,
  # so the expected Fisher X^T W X has a large condition number and the b1
  # posterior is a thin wall -- the structural rehearsal of the eBird tweedie.
  set.seed(202)
  n <- 150L; x <- stats::rnorm(n) * 3
  b_true <- c(1.0, 0.5)
  y <- stats::rpois(n, exp(b_true[1] + b_true[2] * x))
  X <- cbind(1, x)
  stan_code <- paste(
    "data { int<lower=1> n; vector[n] x; array[n] int<lower=0> y; }",
    "parameters { real b0; real b1; }",
    "model { b0 ~ normal(0, 5); b1 ~ normal(0, 5);",
    "        y ~ poisson_log(b0 + b1 * x); }", sep = "\n")
  f <- cmdstanr::write_stan_file(stan_code)
  mod <- cmdstanr::cmdstan_model(f, compile_model_methods = TRUE,
                                 force_recompile = TRUE)
  ref_fit <- mod$sample(data = list(n = n, x = x, y = y), chains = 2,
                        iter_warmup = 1000, iter_sampling = 1000, refresh = 0,
                        show_messages = FALSE, show_exceptions = FALSE)
  ref_draws <- posterior::as_draws_matrix(ref_fit$draws(c("b0", "b1")))
  ref_mean <- colMeans(ref_draws); ref_sd <- apply(ref_draws, 2, stats::sd)
  fit1 <- mod$sample(data = list(n = n, x = x, y = y), chains = 1,
                     iter_warmup = 1, iter_sampling = 1, refresh = 0,
                     show_messages = FALSE, show_exceptions = FALSE)
  eta <- function(b) X %*% b
  tgt <- gdpar_geom_target(
    fit1, dim = 2,
    simulate = function(b) stats::rpois(n, exp(as.numeric(eta(b)))),
    score = function(b, ysim)
      as.numeric(crossprod(X, ysim - exp(as.numeric(eta(b))))))
  fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 8000, seed = 9)
  # the sub-Riemannian geometry built at the posterior mode, walls handled exact.
  metric <- gdpar_geom_metric_subriemannian(tgt, fisher = fisher,
                                            reference = as.numeric(ref_mean))
  expect_gte(metric$n_walls, 1L)                  # at least one stiff wall.
  # the step is matched to the floor (the accessible scale), not the wall: the
  # exact walls impose no limit, which is the whole speed-up over Euclidean HMC.
  eps <- metric$suggested_epsilon
  sub <- gdpar_geom_hmc(tgt, metric = metric, epsilon = eps, L = 25,
                        n_iter = 2000, n_warmup = 800,
                        init = as.numeric(ref_mean), seed = 7)
  expect_true(all(is.finite(sub$draws)))
  expect_lt(sub$n_divergent, 120L)
  sub_mean <- colMeans(sub$draws); sub_sd <- apply(sub$draws, 2, stats::sd)
  # matches NUTS on the original scale (mean and sd of both coefficients).
  expect_lt(max(abs(sub_mean - ref_mean)), 0.05)
  expect_lt(max(abs(sub_sd - ref_sd) / ref_sd), 0.4)
  # at the same step a Euclidean sampler is wrecked by the stiff wall.
  euc <- gdpar_geom_hmc(tgt, epsilon = eps, L = 25, n_iter = 800,
                        n_warmup = 400, init = as.numeric(ref_mean), seed = 7)
  expect_gt(euc$n_divergent, sub$n_divergent + 200L)
})

# ===========================================================================
# RG.4: Finsler / relativistic geometry -- a bounded, non-Gaussian kinetic
# energy coupled to the position-dependent Riemannian metric of RG.3. The remedy
# for heavy tails and directional anisotropy (the G3 target). The velocity
# grad_p K = c M^{-1} p / sqrt(p^T M^{-1} p + m^2 c^2) is capped at the speed c,
# so a large momentum in a heavy tail no longer overshoots; the 0.5 log det M
# normaliser keeps the theta-marginal exact for any c and m. K depends on both
# theta (through M) and p -> the Hamiltonian is non-separable -> a DEDICATED
# generalised implicit leapfrog (Girolami-Calderhead) carried in the integrator
# slot, leaving the default leapfrog bit-identical.
#
# Algebraic layer (ungated, deterministic): kinetic gradients vs finite
# differences; the velocity M-norm is bounded by c; the non-relativistic limit
# (c large -> Gaussian Riemannian kinetic); the dedicated integrator is
# reversible and conserves energy; the momentum sampler reproduces the
# relativistic momentum law; input validation. Statistical layer: ungated (a
# heavy-tailed Student-t recovered through gdpar_geom_hmc) plus gated heavier
# runs (a correlated multivariate Student-t where the Euclidean kinetic
# overshoots, and a real cmdstan-backed correlated multivariate-t matching NUTS).
# ===========================================================================

# --- algebraic layer: the relativistic kinetic energy and its integrator -----

test_that("the relativistic kinetic gradients match finite differences", {
  tgt <- funnel_target(d = 4L, s = 3)              # a position-dependent metric.
  metric <- gdpar_geom_metric_relativistic(
    tgt, fisher = funnel_fisher(4L, 3), dfisher = funnel_dfisher(4L, 3),
    speed = 5, rest_mass = 1)
  kin <- metric$kinetic
  set.seed(11)
  theta <- c(0.3, -0.5, 0.8, -0.2); p <- c(0.7, -0.4, 0.9, 0.2)
  h <- 1e-6
  # grad_p vs central differences of the value in p.
  gp_fd <- vapply(seq_along(p), function(i) {
    e <- numeric(length(p)); e[i] <- h
    (kin$value(theta, p + e) - kin$value(theta, p - e)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(kin$grad_p(theta, p) - gp_fd)), 1e-6)
  # grad_theta vs central differences of the value in theta (through M(theta)).
  gt_fd <- vapply(seq_along(theta), function(i) {
    e <- numeric(length(theta)); e[i] <- h
    (kin$value(theta + e, p) - kin$value(theta - e, p)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(kin$grad_theta(theta, p) - gt_fd)), 1e-5)
})

test_that("the relativistic velocity is M-norm bounded by the speed of light", {
  tgt <- std_normal_target(3L)
  Sigma_inv <- matrix(c(2, 0.5, 0, 0.5, 3, 0.2, 0, 0.2, 1), 3, 3)
  Sigma_inv <- (Sigma_inv + t(Sigma_inv)) / 2
  cspeed <- 4
  metric <- gdpar_geom_metric_relativistic(
    tgt, fisher = function(th) Sigma_inv, speed = cspeed, rest_mass = 1)
  M <- metric$mass(c(0, 0, 0)); kin <- metric$kinetic
  set.seed(7)
  for (scale in c(0.1, 1, 10, 1e3, 1e6)) {            # including enormous momenta.
    p <- scale * stats::rnorm(3)
    v <- kin$grad_p(c(0, 0, 0), p)
    vMv <- as.numeric(crossprod(v, M %*% v))
    expect_lt(sqrt(vMv), cspeed)                        # strictly below c.
  }
})

test_that("the relativistic kinetic reduces to the Gaussian one as c grows", {
  tgt <- std_normal_target(3L)
  Sinv <- diag(c(1, 4, 9))
  metric <- gdpar_geom_metric_relativistic(
    tgt, fisher = function(th) Sinv, speed = 1e4, rest_mass = 1)
  kin <- metric$kinetic; Minv <- metric$inv_mass(c(0, 0, 0))
  set.seed(3); p <- stats::rnorm(3)
  gauss_v <- as.numeric(Minv %*% p)                    # the c -> infinity velocity.
  rel_err <- max(abs(kin$grad_p(c(0, 0, 0), p) - gauss_v)) / max(abs(gauss_v))
  expect_lt(rel_err, 1e-5)
})

test_that("the dedicated relativistic integrator is exactly reversible", {
  tgt <- funnel_target(d = 3L, s = 2.5)
  metric <- gdpar_geom_metric_relativistic(
    tgt, fisher = funnel_fisher(3L, 2.5), dfisher = funnel_dfisher(3L, 2.5),
    speed = 6, rest_mass = 1)
  th0 <- c(0.4, -0.3, 0.5); p0 <- c(-0.6, 0.2, 0.35)
  fwd <- metric$integrator(th0, p0, tgt, 0.25, 8L)
  rev <- metric$integrator(fwd$theta, -fwd$p, tgt, 0.25, 8L)
  expect_true(fwd$converged && rev$converged)
  expect_lt(max(abs(rev$theta - th0)), 1e-8)           # position returns...
  expect_lt(max(abs(-rev$p - p0)), 1e-8)               # ...and so does the momentum.
})

test_that("the relativistic integrator conserves energy without drift", {
  # A constant dense mass + Gaussian potential: a relativistic-Gaussian system on
  # which the implicit leapfrog conserves H to second order (bounded, no drift).
  tgt <- gdpar_geom_target(
    log_prob = function(th) -0.5 * (th[1]^2 + 6 * th[2]^2),
    grad_log_prob = function(th) -c(th[1], 6 * th[2]), dim = 2)
  metric <- gdpar_geom_metric_relativistic(
    tgt, fisher = function(th) diag(c(1, 6)), speed = 3, rest_mass = 1)
  kin <- metric$kinetic
  H <- function(th, pp) -tgt$log_prob(th) + kin$value(th, pp)
  th0 <- c(1.1, 0.4); p0 <- kin$draw_momentum(th0)
  H0 <- H(th0, p0)
  drift <- vapply(1:25, function(L) {
    st <- metric$integrator(th0, p0, tgt, 0.1, L)
    abs(H(st$theta, st$p) - H0)
  }, numeric(1))
  expect_lt(max(drift), 5e-3)                           # bounded, no secular growth.
})

test_that("the relativistic momentum sampler reproduces the relativistic law", {
  tgt <- std_normal_target(3L)
  cspeed <- 2.5; rmass <- 1
  metric <- gdpar_geom_metric_relativistic(
    tgt, fisher = function(th) diag(3), speed = cspeed, rest_mass = rmass)
  kin <- metric$kinetic
  # theoretical radial second moment of r ~ r^(d-1) exp(-c sqrt(r^2 + m^2 c^2)).
  d <- 3L; mc2 <- (rmass * cspeed)^2
  rr <- seq(1e-6, 60, length.out = 200000L)
  lw <- (d - 1) * log(rr) - cspeed * sqrt(rr^2 + mc2)
  w <- exp(lw - max(lw))
  Er2 <- sum(rr^2 * w) / sum(w)
  set.seed(99)
  P <- t(vapply(1:6000, function(i) kin$draw_momentum(c(0, 0, 0)), numeric(3)))
  expect_lt(abs(mean(rowSums(P^2)) - Er2) / Er2, 0.05)  # radius^2 moment matches.
  Cov <- stats::cov(P)                                  # isotropic (mass = I).
  expect_lt(max(abs(Cov[upper.tri(Cov)])) / (Er2 / d), 0.12)
  expect_lt(max(abs(diag(Cov) - Er2 / d)) / (Er2 / d), 0.15)
})

test_that("gdpar_geom_metric_relativistic exposes its structure and validates input", {
  tgt <- std_normal_target(2L)
  metric <- gdpar_geom_metric_relativistic(tgt, fisher = function(th) diag(2),
                                           speed = 7, rest_mass = 2)
  expect_identical(metric$metric_kind, "relativistic")
  expect_true(metric$position_dependent)
  expect_equal(metric$speed, 7); expect_equal(metric$rest_mass, 2)
  expect_true(is.function(metric$integrator) && is.list(metric$kinetic))
  expect_error(gdpar_geom_metric_relativistic(tgt, fisher = function(th) diag(2),
                                              speed = 0),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_metric_relativistic(tgt, fisher = function(th) diag(2),
                                              rest_mass = -1),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_metric_relativistic(tgt, curvature = "fisher"),
               class = "gdpar_input_error")          # fisher required.
})

# --- statistical layer (ungated, light): a heavy-tailed Student-t ------------

test_that("the relativistic sampler recovers a heavy-tailed Student-t", {
  nu <- 4; d <- 2L
  tgt <- gdpar_geom_target(
    log_prob = function(th) -((nu + 1) / 2) * sum(log1p(th^2 / nu)),
    grad_log_prob = function(th) -(nu + 1) * th / (nu + th^2), dim = d)
  fisher <- function(th) diag((nu + 1) / (nu + 3), d)   # expected Fisher of the t.
  metric <- gdpar_geom_metric_relativistic(tgt, fisher = fisher, speed = 4,
                                           rest_mass = 1)
  fit <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.5, L = 12,
                        n_iter = 3000, n_warmup = 800, seed = 5)
  expect_true(all(is.finite(fit$draws)))
  expect_gt(fit$accept_rate, 0.6)
  expect_lt(fit$n_divergent, 30L)
  # the marginal variance of a Student-t is nu / (nu - 2) = 2 for nu = 4.
  expect_lt(abs(stats::var(fit$draws[, 1]) - nu / (nu - 2)) / (nu / (nu - 2)), 0.4)
})

# --- statistical layer (gated): correlated multivariate-t (tails + anisotropy)

test_that("[gated] relativistic samples a correlated multivariate-t where Euclidean overshoots", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the multivariate-t run.")
  }
  d <- 4L; nu <- 5
  R <- matrix(0.6, d, d); diag(R) <- 1                  # directional correlation.
  s <- c(1, 2, 0.5, 1.5)                                # anisotropic scales.
  Sigma <- diag(s) %*% R %*% diag(s); Sinv <- solve(Sigma)
  lp <- function(th) -((nu + d) / 2) *
    log1p(as.numeric(crossprod(th, Sinv %*% th)) / nu)
  gr <- function(th) {
    q <- as.numeric(Sinv %*% th)
    -(nu + d) * q / (nu + as.numeric(crossprod(th, q)))
  }
  tgt <- gdpar_geom_target(log_prob = lp, grad_log_prob = gr, dim = d)
  fisher <- function(th) ((nu + d) / (nu + d + 2)) * Sinv  # expected Fisher (const).
  metric <- gdpar_geom_metric_relativistic(tgt, fisher = fisher, speed = 6,
                                           rest_mass = 1)
  rel <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.6, L = 16,
                        n_iter = 6000, n_warmup = 1500, seed = 2)
  expect_true(all(is.finite(rel$draws)))
  expect_lt(rel$n_divergent, 60L)
  expect_gt(rel$accept_rate, 0.6)
  # the covariance of a multivariate-t is nu/(nu-2) Sigma; recover it.
  Cov_true <- nu / (nu - 2) * Sigma
  Cov_hat <- stats::cov(rel$draws)
  expect_lt(max(abs(Cov_hat - Cov_true)) / max(abs(Cov_true)), 0.35)
  # at the same step a Euclidean (identity-mass Gaussian) kinetic overshoots the
  # anisotropic heavy tails: its energy errors are far larger, so its acceptance
  # collapses (the bounded velocity is what keeps the relativistic proposal
  # informative). The clean signal here is acceptance, not divergence (the
  # overshoots are large-energy rejections, not |dH| > 1000 divergences).
  euc <- gdpar_geom_hmc(tgt, epsilon = 0.6, L = 16, n_iter = 6000,
                        n_warmup = 1500, seed = 2)
  expect_gt(rel$accept_rate, euc$accept_rate + 0.2)
})

# --- statistical layer (gated, real cmdstan): correlated multivariate-t ------

test_that("[gated] relativistic-Riemannian matches NUTS on a cmdstan multivariate-t", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the cmdstan multivariate-t demo.")
  }
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")
  d <- 3L; nu <- 4
  R <- matrix(0.5, d, d); diag(R) <- 1
  s <- c(1, 1.8, 0.6); Sigma <- diag(s) %*% R %*% diag(s); Sinv <- solve(Sigma)
  stan_code <- paste(
    "data { int<lower=1> d; real<lower=1> nu; vector[d] mu; matrix[d, d] Sigma; }",
    "parameters { vector[d] theta; }",
    "model { theta ~ multi_student_t(nu, mu, Sigma); }", sep = "\n")
  f <- cmdstanr::write_stan_file(stan_code)
  mod <- cmdstanr::cmdstan_model(f, compile_model_methods = TRUE,
                                 force_recompile = TRUE)
  sdata <- list(d = d, nu = nu, mu = rep(0, d), Sigma = Sigma)
  ref_fit <- mod$sample(data = sdata, chains = 2, iter_warmup = 1000,
                        iter_sampling = 1500, refresh = 0,
                        show_messages = FALSE, show_exceptions = FALSE)
  ref_draws <- posterior::as_draws_matrix(
    ref_fit$draws(paste0("theta[", 1:d, "]")))
  ref_mean <- colMeans(ref_draws); ref_cov <- stats::cov(ref_draws)
  fit1 <- mod$sample(data = sdata, chains = 1, iter_warmup = 1, iter_sampling = 1,
                     refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  tgt <- gdpar_geom_target(fit1, dim = d)
  # The bounded relativistic kinetic on the natural (constant) expected Fisher of
  # the multivariate-t, (nu+d)/(nu+d+2) Sigma^{-1}, driving the REAL cmdstan
  # log-density / gradient: a real heavy-tailed posterior sampled exactly where a
  # Gaussian kinetic overshoots the tails. (The SoftAbs mass is the funnel remedy,
  # mismatched to an elliptical heavy-tailed target; the coupled
  # position-dependent integrator is validated on the funnel below and above.)
  fisher <- function(th) ((nu + d) / (nu + d + 2)) * Sinv
  metric <- gdpar_geom_metric_relativistic(tgt, fisher = fisher, speed = 5,
                                           rest_mass = 1)
  rel <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.5, L = 16,
                        n_iter = 2500, n_warmup = 1000, init = as.numeric(ref_mean),
                        seed = 4)
  expect_true(all(is.finite(rel$draws)))
  expect_lt(rel$n_divergent, 60L)
  expect_gt(rel$accept_rate, 0.7)
  rel_mean <- colMeans(rel$draws); rel_cov <- stats::cov(rel$draws)
  expect_lt(max(abs(rel_mean - ref_mean)), 0.2)           # matches NUTS location...
  expect_lt(max(abs(rel_cov - ref_cov)) / max(abs(ref_cov)), 0.4)  # ...and spread.
})

# --- statistical layer (gated): the coupled, position-dependent case ---------
# Neal's funnel exercises the non-separable integrator end to end: a
# position-dependent Fisher mass (variable curvature) under the bounded
# relativistic kinetic. The relativistic-Riemannian sampler resolves the neck
# where the Euclidean kinetic stalls -- the coupling the user chose (funnel and
# heavy tails together), here on the funnel half.

test_that("[gated] coupled relativistic-Riemannian resolves Neal's funnel", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the coupled funnel run.")
  }
  d <- 10L; s <- 3                                  # true sd(v) = 3.
  tgt <- funnel_target(d, s)
  metric <- gdpar_geom_metric_relativistic(
    tgt, curvature = "fisher", fisher = funnel_fisher(d, s),
    dfisher = funnel_dfisher(d, s), speed = 10, rest_mass = 1)
  rel <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.3, L = 12,
                        n_iter = 800, n_warmup = 400, seed = 42)
  euc <- gdpar_geom_hmc(tgt, epsilon = 0.1, L = 20, n_iter = 800,
                        n_warmup = 400, seed = 42)
  expect_lt(rel$n_divergent, 20L)
  expect_gt(rel$accept_rate, 0.7)
  expect_gt(rel$ebfmi, 0.5)                          # the metric resolves the neck.
  expect_lt(euc$ebfmi, 0.2)                          # the Euclidean metric does not.
  expect_gt(rel$ebfmi, 3 * euc$ebfmi)
  expect_lt(min(rel$draws[, 1]), min(euc$draws[, 1]) - 1)   # reaches deeper.
  expect_lt(abs(stats::sd(rel$draws[, 1]) - s), 0.4)        # recovers Var(v).
})

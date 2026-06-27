# Deep, gated reproduction of the heavy (cmdstan-backed) runs narrated in the
# Block RG vignette `vop08_geometric_robustness.Rmd`. These compile Stan models
# and sample, so they take several minutes and are NEVER run during a normal
# package check. Run them on purpose:
#
#   Sys.setenv(GDPAR_RUN_GEOMETRY_PILOTS = "1")
#   source(system.file("scripts", "geometry_pilots_deep.R", package = "gdpar"))
#
# Everything here is opt-in and read-only with respect to the package: it never
# touches gdpar()'s default fit path. The vignette evaluates only the cheap,
# closure-based counterparts of these runs; this script is the full,
# falsifiable end-to-end reproduction.

if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
  message("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to run the deep geometry pilots; skipping.")
} else if (!requireNamespace("cmdstanr", quietly = TRUE) ||
           !requireNamespace("posterior", quietly = TRUE)) {
  message("cmdstanr and posterior are required for the deep geometry pilots; skipping.")
} else {
  library(gdpar)
  suite <- gdpar_geometry_suite()
  rule <- function(s) cat("\n========== ", s, " ==========\n", sep = "")

  ## (1) The diagnostic on a compiled diana ---------------------------------
  rule("RG.1  diagnostic on a compiled cmdstan diana (G4 near-deterministic)")
  diag <- gdpar_geometry_diagnostic(suite$G4_quasi_deterministic, n_grid = 3)
  print(diag)
  cat("pathology:", diag$pathology, " culprit:",
      paste(diag$culprit$name, collapse = ", "), "\n")

  ## (2) Each ladder level on a real cmdstan model --------------------------
  # A compiled Neal funnel (position-dependent curvature -> Riemannian level).
  rule("RG.2-RG.3  Euclidean vs Riemannian (SoftAbs) on a compiled funnel")
  funnel_code <- "
data { int<lower=1> K; real s; }
parameters { real v; vector[K] x; }
model { v ~ normal(0, s); x ~ normal(0, exp(v / 2)); }
"
  fm <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(funnel_code),
                                compile_model_methods = TRUE)
  ff <- fm$sample(data = list(K = 9L, s = 3), chains = 1, iter_warmup = 1,
                  iter_sampling = 1, refresh = 0, show_messages = FALSE,
                  show_exceptions = FALSE)
  ff$init_model_methods(verbose = FALSE, hessian = TRUE)
  d_funnel <- ncol(posterior::as_draws_matrix(ff$unconstrain_draws()))
  tgt_f <- gdpar_geom_target(object = ff, dim = d_funnel)
  eucl <- gdpar_geom_hmc(tgt_f, gdpar_geom_metric_euclidean(dim = d_funnel),
                         epsilon = 0.1, L = 20, n_iter = 800, n_warmup = 400,
                         seed = 1)
  rie <- gdpar_geom_hmc(tgt_f, gdpar_geom_metric_riemannian(tgt_f, "softabs"),
                        epsilon = 0.3, L = 12, n_iter = 800, n_warmup = 400,
                        seed = 1)
  cat(sprintf("Euclidean  E-BFMI %.3f accept %.3f\n", eucl$ebfmi, eucl$accept_rate))
  cat(sprintf("Riemannian E-BFMI %.3f accept %.3f\n", rie$ebfmi, rie$accept_rate))

  ## (3) The full orchestrator: resolution and a certified limit ------------
  rule("RG.5  orchestrator -- resolution on the isotropic control")
  b <- gdpar_geom_orchestrate_budget()
  b$tune_epsilon <- FALSE; b$probe_iter <- 80L; b$full_iter <- 100L
  b$full_warmup <- 100L
  res <- gdpar_geom_orchestrate(suite$G0_isotropic, n_grid = 1, budget = b,
                                pilot_warmup = 100L, pilot_sampling = 100L,
                                verbose = FALSE)
  print(res)

  rule("RG.5  orchestrator -- certified limit when the Fisher is withheld")
  # A near-deterministic target whose remedy (sub-Riemannian) needs a Fisher:
  # without one the controller certifies the limit and prescribes supplying it.
  res_cert <- gdpar_geom_orchestrate(suite$G4_quasi_deterministic, n_grid = 1,
                                     fisher = NULL, budget = b, verbose = FALSE)
  print(res_cert)
  if (!is.null(res_cert$certificate)) print(res_cert$certificate)

  ## (4) The bridge over a real fit + the one-call gdpar_geom_fit ------------
  rule("RG.6 (ii)  bridge over a fitted gdpar model + gdpar_geom_fit")
  set.seed(1)
  n <- 80
  x <- rnorm(n); z <- rnorm(n)
  y <- rnorm(n, 0.5 + 0.8 * (x - mean(x)), exp(-0.2 + 0.4 * (z - mean(z))))
  dat <- data.frame(y = y, x = x, z = z)
  fit <- gdpar(gdpar_bf(y ~ a(x), sigma ~ a(z)), data = dat,
               family = gdpar_family("gaussian"), chains = 1,
               iter_warmup = 200, iter_sampling = 200, refresh = 0,
               seed = 313, skip_id_check = TRUE, verbose = FALSE)
  bridge <- gdpar_geom_bridge(fit)
  print(bridge)
  b2 <- gdpar_geom_orchestrate_budget()
  b2$tune_epsilon <- FALSE; b2$max_rounds <- 2L
  b2$probe_iter <- 80L; b2$full_iter <- 100L; b2$full_warmup <- 100L
  res_bridge <- gdpar_geom_orchestrate(bridge$target, bridge$geom_target,
                                       reference = bridge$reference,
                                       budget = b2, n_grid = 1, verbose = FALSE)
  print(res_bridge)

  res_fit <- gdpar_geom_fit(gdpar_bf(y ~ a(x), sigma ~ a(z)), data = dat,
                            family = gdpar_family("gaussian"),
                            skip_id_check = TRUE, budget = b2, n_grid = 1,
                            verbose = FALSE)
  print(res_fit)
  rule("DONE -- deep geometry pilots complete")
}

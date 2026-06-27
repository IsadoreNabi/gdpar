# Block RG, RG.5: geometry-adaptive sampling orchestrator.
#
# Default (ungated) tests are pure R (no cmdstan): the tunable factories, the
# ladder/selection/gate/score/prescription/certificate logic and the armour
# (monotone ratchet helpers, deterministic seeding, input validation), exercised
# with synthetic diagnostics and fits. The end-to-end layer -- orchestrating the
# real suite targets and a compiled cmdstan model through the diagnostic pilots
# -- is gated by GDPAR_RUN_GEOMETRY_PILOTS, the same gate as the other RG pilots.

# --- helpers ---------------------------------------------------------------

# A synthetic gdpar_geometry_diagnostic-like list carrying exactly the fields the
# orchestrator's selection / prescription helpers read.
mock_diag <- function(pathology, confidence, condition_number, heavy_kurtosis,
                      ebfmi_min, grows_with_n, slope = 1, tractability = "tractable") {
  list(
    pathology = pathology, confidence = confidence,
    signals = data.frame(condition_number = condition_number,
                         heavy_kurtosis = heavy_kurtosis, ebfmi_min = ebfmi_min),
    difficulty_curve = list(slope = slope, grows_with_n = grows_with_n),
    culprit = data.frame(parameter = "v", mechanism = "neck", score = 1,
                         stringsAsFactors = FALSE),
    cost = list(tractability = tractability, seconds_per_1000_draws = 1))
}

# --- tunable factories -----------------------------------------------------

test_that("the budget and criteria factories expose tunable data", {
  b <- gdpar_geom_orchestrate_budget()
  cr <- gdpar_geom_orchestrate_criteria()
  expect_true(is.list(b) && b$max_rounds >= 1L && b$max_levels >= 1L)
  expect_true(all(c("probe_iter", "full_iter", "hysteresis", "stall_limit",
                    "max_fits") %in% names(b)))
  expect_true(is.list(cr) && cr$accept_low < cr$accept_high)
  expect_true(cr$ebfmi_low > 0 && cr$divergent_rate_high > 0)
})

# --- ladder and monotone climb --------------------------------------------

test_that("the pathology-to-level map and the monotone climb behave", {
  expect_identical(gdpar:::.gdpar_orch_pathology_to_key("isotropic"),
                   "euclidean_diagonal")
  expect_identical(gdpar:::.gdpar_orch_pathology_to_key("funnel"), "riemannian")
  expect_identical(gdpar:::.gdpar_orch_pathology_to_key("heavy_tails"),
                   "relativistic")
  expect_identical(gdpar:::.gdpar_orch_pathology_to_key("quasi_deterministic"),
                   "sub_riemannian")
  # Out-of-scope remedies map to NA (handled by a short-circuit certificate).
  expect_true(is.na(gdpar:::.gdpar_orch_pathology_to_key("multimodal")))
  expect_true(is.na(gdpar:::.gdpar_orch_pathology_to_key("flat_direction")))
  # The ladder climbs strictly upward and terminates.
  expect_identical(gdpar:::.gdpar_orch_next_above(0L), "euclidean_dense")
  expect_identical(gdpar:::.gdpar_orch_next_above(1L), "riemannian")
  expect_identical(gdpar:::.gdpar_orch_next_above(3L), "relativistic")
  expect_identical(gdpar:::.gdpar_orch_next_above(4L), "sub_riemannian")
  expect_true(is.na(gdpar:::.gdpar_orch_next_above(5L)))
})

# --- deterministic seeding -------------------------------------------------

test_that("the seeding schedule is deterministic and slot-sensitive", {
  s <- gdpar:::.gdpar_orch_seed(20260603L, 2L, 1L)
  expect_identical(s, gdpar:::.gdpar_orch_seed(20260603L, 2L, 1L))
  expect_true(s >= 1L)
  expect_false(identical(gdpar:::.gdpar_orch_seed(1L, 1L, 1L),
                         gdpar:::.gdpar_orch_seed(1L, 1L, 2L)))
  expect_false(identical(gdpar:::.gdpar_orch_seed(1L, 1L, 1L),
                         gdpar:::.gdpar_orch_seed(1L, 2L, 1L)))
})

# --- (a) combined selection -------------------------------------------------

test_that("a clear, confident pathology uses the discrete level directly", {
  th <- gdpar_geometry_thresholds()
  d_fun <- mock_diag("funnel", 0.9, condition_number = 15, heavy_kurtosis = 4,
                     ebfmi_min = 0.1, grows_with_n = FALSE)
  sel <- gdpar:::.gdpar_orch_select_entry(d_fun, NULL, NULL, th)
  expect_identical(sel$key, "riemannian")
  expect_false(sel$conflict)
  d_q <- mock_diag("quasi_deterministic", 0.9, condition_number = 50,
                   heavy_kurtosis = 0, ebfmi_min = 0.9, grows_with_n = TRUE,
                   slope = 1.2)
  expect_identical(gdpar:::.gdpar_orch_select_entry(d_q, NULL, NULL, th)$key,
                   "sub_riemannian")
})

test_that("a discrete/proximity conflict starts at the lower (conservative) level", {
  th <- gdpar_geometry_thresholds()
  # Discrete says heavy_tails (level 4) but the depressed E-BFMI is the funnel
  # signature (level 3); the controller starts conservatively at riemannian (3).
  d_c <- mock_diag("heavy_tails", 0.4, condition_number = 5, heavy_kurtosis = 4,
                   ebfmi_min = 0.05, grows_with_n = FALSE)
  sel <- gdpar:::.gdpar_orch_select_entry(d_c, NULL, NULL, th)
  expect_true(sel$conflict)
  expect_identical(sel$key, "riemannian")
})

test_that("an out-of-scope pathology yields an NA entry key", {
  th <- gdpar_geometry_thresholds()
  d_m <- mock_diag("multimodal", 0.9, condition_number = 2, heavy_kurtosis = 0,
                   ebfmi_min = 0.9, grows_with_n = FALSE)
  expect_true(is.na(gdpar:::.gdpar_orch_select_entry(d_m, NULL, NULL, th)$key))
})

test_that("user overrides take precedence over the combined map", {
  th <- gdpar_geometry_thresholds()
  d_fun <- mock_diag("funnel", 0.9, condition_number = 15, heavy_kurtosis = 4,
                     ebfmi_min = 0.1, grows_with_n = FALSE)
  s_map <- gdpar:::.gdpar_orch_select_entry(d_fun, list(funnel = "relativistic"),
                                            NULL, th)
  expect_identical(s_map$key, "relativistic")
  expect_identical(s_map$source, "user_level_map")
  s_ent <- gdpar:::.gdpar_orch_select_entry(d_fun, NULL, "euclidean_dense", th)
  expect_identical(s_ent$key, "euclidean_dense")
  expect_error(gdpar:::.gdpar_orch_select_entry(d_fun, NULL, "nonsense", th),
               class = "gdpar_input_error")
})

# --- multi-signal success gate + score -------------------------------------

test_that("the success gate is a conjunction with hysteresis", {
  cr <- gdpar_geom_orchestrate_criteria()
  good <- list(accept_rate = 0.85, n_divergent = 0L, ebfmi = 0.95)
  bad <- list(accept_rate = 0.2, n_divergent = 200L, ebfmi = 0.02)
  expect_true(gdpar:::.gdpar_orch_success_gate(good, 300L, cr, 0)$pass)
  expect_false(gdpar:::.gdpar_orch_success_gate(bad, 300L, cr, 0.1)$pass)
  # A marginal fit passes the bare probe gate but fails the tightened full gate.
  marg <- list(accept_rate = 0.85, n_divergent = 0L, ebfmi = 0.31)
  expect_true(gdpar:::.gdpar_orch_success_gate(marg, 300L, cr, 0)$pass)
  expect_false(gdpar:::.gdpar_orch_success_gate(marg, 300L, cr, 0.1)$pass)
  # A high-acceptance / zero-mixing fit must NOT pass (the rhat=Inf lesson).
  stuck <- list(accept_rate = 1, n_divergent = 0L, ebfmi = 0.01)
  expect_false(gdpar:::.gdpar_orch_success_gate(stuck, 300L, cr, 0)$pass)
})

test_that("the health score orders a good fit above a bad one", {
  good <- list(accept_rate = 0.85, n_divergent = 0L, ebfmi = 0.95)
  bad <- list(accept_rate = 0.2, n_divergent = 200L, ebfmi = 0.02)
  expect_gt(gdpar:::.gdpar_orch_score(good, 300L),
            gdpar:::.gdpar_orch_score(bad, 300L))
  expect_identical(gdpar:::.gdpar_orch_score(NULL, 300L), -Inf)
})

# --- (c) prescription + certificate ----------------------------------------

test_that("the prescription is evidence-driven and tagged as conjectured", {
  b <- gdpar_geom_orchestrate_budget()
  d_q <- mock_diag("quasi_deterministic", 0.9, condition_number = 50,
                   heavy_kurtosis = 0, ebfmi_min = 0.9, grows_with_n = TRUE,
                   slope = 1.2)
  presc <- gdpar:::.gdpar_orch_prescription(d_q, list(), FALSE, b, 0.5,
                                            character(0))
  acts <- vapply(presc, function(p) p$action, character(1))
  expect_true(any(grepl("expected Fisher", acts)))
  expect_true(all(vapply(presc, function(p) identical(p$status, "conjectured"),
                         logical(1))))
  expect_true(all(vapply(presc, function(p) nzchar(p$falsifiable_test),
                         logical(1))))
  d_m <- mock_diag("multimodal", 0.9, condition_number = 2, heavy_kurtosis = 0,
                   ebfmi_min = 0.9, grows_with_n = FALSE)
  expect_true(any(grepl("tempering",
    vapply(gdpar:::.gdpar_orch_prescription(d_m, list(), TRUE, b, -Inf,
                                            character(0)),
           function(p) p$action, character(1)))))
})

test_that("the certificate separates the three rigour layers", {
  b <- gdpar_geom_orchestrate_budget()
  d_q <- mock_diag("quasi_deterministic", 0.9, condition_number = 50,
                   heavy_kurtosis = 0, ebfmi_min = 0.9, grows_with_n = TRUE)
  presc <- gdpar:::.gdpar_orch_prescription(d_q, list(), FALSE, b, 0.5,
                                            character(0))
  cert <- gdpar:::.gdpar_orch_certificate("intractable_within_budget", d_q,
            list(), d_q$culprit, presc, list(seed = 1),
            list(fits = 3L, seconds = 10))
  expect_s3_class(cert, "gdpar_geom_certificate")
  expect_true(!is.null(cert$algebraic) && !is.null(cert$statistical) &&
              !is.null(cert$numerical))
  expect_identical(cert$algebraic$grows_with_n, TRUE)
  expect_output(print(cert), "gdpar_geom_certificate")
})

# --- input validation ------------------------------------------------------

test_that("gdpar_geom_orchestrate validates its arguments before any sampling", {
  suite <- gdpar_geometry_suite()
  expect_error(gdpar_geom_orchestrate(suite$G0_isotropic, fisher = 5,
                                      verbose = FALSE),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_orchestrate(suite$G0_isotropic,
                                      reference = c(1, 2, 3), verbose = FALSE),
               class = "gdpar_input_error")
  expect_error(gdpar_geom_orchestrate(suite$G0_isotropic, speed = -1,
                                      verbose = FALSE),
               class = "gdpar_input_error")
})

# --- gated end-to-end layer (cmdstan pilots + R-native sampling) ------------

# A small, fast budget for the gated routing tests.
orch_budget_fast <- function() {
  b <- gdpar_geom_orchestrate_budget()
  b$probe_warmup <- 80L; b$probe_iter <- 80L
  b$full_warmup <- 120L; b$full_iter <- 120L
  b$tune_iter <- 40L; b$L <- 18L; b$max_rounds <- 4L
  b
}

test_that("the orchestrator routes the easy and contracting targets correctly", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the orchestration routing.")
  }
  suite <- gdpar_geometry_suite()
  b <- orch_budget_fast()
  # Isotropic control resolves at the cheapest Euclidean level.
  r0 <- gdpar_geom_orchestrate(suite$G0_isotropic, n_grid = c(1, 4), budget = b,
                               pilot_warmup = 80L, pilot_sampling = 80L,
                               verbose = FALSE)
  expect_identical(r0$status, "resolved")
  expect_identical(r0$level, "euclidean_diagonal")
  # Quasi-deterministic canyon, with the expected Fisher, resolves at the
  # sub-Riemannian level.
  n_eff <- suite$G4_quasi_deterministic$default_n
  dd <- suite$G4_quasi_deterministic$default_difficulty
  r4 <- gdpar_geom_orchestrate(suite$G4_quasi_deterministic,
          n_grid = c(10, 100, 1000),
          fisher = function(theta) diag(c(1, rep(n_eff, dd - 1))),
          budget = b, pilot_warmup = 80L, pilot_sampling = 80L, verbose = FALSE)
  expect_identical(r4$status, "resolved")
  expect_identical(r4$level, "sub_riemannian")
})

test_that("the orchestrator certifies a limit with an actionable Fisher prescription", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the orchestration routing.")
  }
  suite <- gdpar_geometry_suite()
  b <- orch_budget_fast()
  # The same canyon WITHOUT the Fisher: the sub-Riemannian level is unbuildable,
  # so a certified limit is emitted, and its conjectured prescription names the
  # smallest fix -- supply the expected Fisher (which, supplied, resolves it).
  r <- gdpar_geom_orchestrate(suite$G4_quasi_deterministic,
         n_grid = c(10, 100, 1000), budget = b, pilot_warmup = 80L,
         pilot_sampling = 80L, verbose = FALSE)
  expect_identical(r$status, "certified_limit")
  acts <- vapply(r$certificate$prescription, function(p) p$action, character(1))
  expect_true(any(grepl("expected Fisher", acts)))
  expect_true(all(vapply(r$certificate$prescription,
                         function(p) identical(p$status, "conjectured"),
                         logical(1))))
})

test_that("the orchestrator routes a multimodal target out of scope", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the orchestration routing.")
  }
  suite <- gdpar_geometry_suite()
  r <- gdpar_geom_orchestrate(suite$G5_multimodal, n_grid = c(1, 4),
         budget = orch_budget_fast(), pilot_warmup = 80L, pilot_sampling = 80L,
         verbose = FALSE)
  expect_identical(r$status, "out_of_scope")
  acts <- vapply(r$certificate$prescription, function(p) p$action, character(1))
  expect_true(any(grepl("tempering", acts)))
  # No spurious geometry advice for an out-of-scope diagnosis.
  expect_false(any(grepl("expected Fisher", acts)))
})

test_that("the orchestrator drives a real compiled cmdstan model (the bridge)", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the cmdstan bridge.")
  }
  skip_if_not_installed("cmdstanr")
  fcode <- paste("data { int<lower=1> d; real<lower=0> scale_v; }",
    "parameters { real v; vector[d - 1] x; }",
    "model { v ~ normal(0, scale_v); x ~ normal(0, exp(v / 2)); }", sep = "\n")
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(fcode),
                                 compile_model_methods = TRUE,
                                 force_recompile = TRUE)
  tgt <- list(model = mod, data = list(d = 3L, scale_v = 2), dim = 3L)
  ff <- function(theta) { v <- theta[1]; diag(c(0.5 * 2 + 1 / 4, rep(exp(-v), 2))) }
  r <- gdpar_geom_orchestrate(tgt, n_grid = c(1L), fisher = ff,
         entry_level = "riemannian", budget = orch_budget_fast(),
         pilot_warmup = 80L, pilot_sampling = 80L, verbose = FALSE)
  expect_identical(r$status, "resolved")
  expect_identical(r$level, "riemannian")
})

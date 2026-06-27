# Block 9.1 -- Internal synthetic adversarial re-validation.
# Layer 1 (Decision C = C.iv, Session B9.7, 2026-05-28): STRUCTURAL
# adversarial fuzz over the canonical architecture. Every test in this
# file deliberately violates one canonical invariant and asserts that
# the package refuses the degenerate input with a canonical abort (or a
# canonical "passed = FALSE" verdict from the identifiability machinery)
# rather than silently producing a fit.
#
# These tests are FAST: every abort happens pre-fit (input validation,
# identifiability algebra, prior-structure validation), so the file runs
# by default in the suite with no Stan compilation.
#
# Naming note: the Charter refers to a "gdpar_c7_violation" class; the
# canonical abort surface actually raises "gdpar_input_error" or
# "gdpar_identifiability_error", and the C4/C4-bis identifiability checks
# return a report with passed = FALSE rather than aborting. The tests
# below assert the *actual* canonical behaviour.

# ---- Section A: C7 / grouping aliased (collinear design) ----------------

test_that("ADV-A1: perfectly aliased additive covariate fails C4 identifiability", {
  df <- make_gaussian_data(n = 200, seed = 42)
  df$x_dup <- df$x1
  spec <- amm_spec(a = ~ x1 + x_dup)
  rep <- gdpar_check_identifiability(spec, df, formula_rhs = ~ x1 + x_dup)
  expect_false(rep$passed)
  expect_false(is.null(rep$collinear_directions))
  expect_gt(rep$condition_number, 1e6)
})

test_that("ADV-A2: scaled-replica covariate (alias up to a constant) still fails C4", {
  df <- make_gaussian_data(n = 200, seed = 7)
  df$x_scaled <- 1000 * df$x1
  spec <- amm_spec(a = ~ x1 + x_scaled)
  rep <- gdpar_check_identifiability(spec, df, formula_rhs = ~ x1 + x_scaled)
  # Column normalization cannot rescue a deterministic linear dependence.
  expect_false(rep$passed)
  expect_false(is.null(rep$collinear_directions))
})

test_that("ADV-A3: degree-1 W aliased with additive block fails per-coord C4-bis", {
  set.seed(42)
  n <- 300
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "polynomial", degree = 1),
    x_vars = c("x1", "x2")
  )
  rep <- gdpar_check_identifiability(
    spec, df, theta_ref_init = c(0.5, -0.3),
    formula_rhs = ~ x1 + x2, rigor = "full"
  )
  expect_false(rep$passed)
  for (pk in rep$c4_bis$per_k) {
    expect_false(pk$passed)
    expect_false(is.null(pk$collinear_directions))
  }
})

# ---- Section B: D-ID per-slot adversarial -------------------------------

test_that("ADV-B1: duplicate prior_canonical_kind across K>=2 slots fails D-ID", {
  fam <- gdpar_family("gaussian")
  # Promote the dispersion slot to a per-group distributional parameter
  # and collide its canonical kind with the location slot's kind.
  fam$param_specs[[2L]]$scope <- "per_group"
  fam$param_specs[[2L]]$prior_canonical_kind <- "mu"
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  rep <- gdpar_check_identifiability(spec, df, family = fam, rigor = "full")
  expect_false(rep$did_pre_fit$symbolic_separability$passed)
  expect_false(rep$did_pre_fit$passed)
  expect_false(rep$passed)
})

test_that("ADV-B2: clean K>=2 family with distinct canonical kinds passes D-ID", {
  fam <- gdpar_family("gaussian")
  fam$param_specs[[2L]]$scope <- "per_group"
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  rep <- gdpar_check_identifiability(spec, df, family = fam, rigor = "full")
  expect_equal(rep$did_pre_fit$K, 2L)
  expect_true(rep$did_pre_fit$symbolic_separability$passed)
  expect_true(rep$did_pre_fit$passed)
})

# ---- Section C: param_specs perturbation (corrupted family object) ------

test_that("ADV-C1: gdpar_family rejects an unknown family name", {
  expect_error(
    gdpar_family("not_a_real_family"),
    regexp = "should be one of"
  )
})

test_that("ADV-C2: the 12 canonical families expose the expected stan_id contract", {
  expected <- c(gaussian = 1L, poisson = 2L, neg_binomial_2 = 3L,
                bernoulli = 4L, beta = 5L, gamma = 6L, student_t = 8L,
                tweedie = 9L, zip = 10L, zinb = 11L, hurdle_poisson = 12L,
                hurdle_neg_binomial_2 = 13L)
  for (nm in names(expected)) {
    fam <- gdpar_family(nm)
    expect_identical(as.integer(fam$stan_id), expected[[nm]],
                     info = sprintf("stan_id drift for %s", nm))
    expect_true(length(fam$param_specs) >= 1L)
    # Every param spec must carry a canonical-kind tag for D-ID.
    for (ps in fam$param_specs) {
      expect_true(is.character(ps$prior_canonical_kind))
      expect_true(nzchar(ps$prior_canonical_kind))
    }
  }
})

# ---- Section D: degenerate priors ---------------------------------------

test_that("ADV-D1: structurally degenerate prior args abort with gdpar_input_error", {
  expect_error(gdpar_prior(sigma_y = ""), class = "gdpar_input_error")
  expect_error(gdpar_prior(sigma_y = c("a", "b")), class = "gdpar_input_error")
  expect_error(gdpar_prior(theta_ref = 0), class = "gdpar_input_error")
})

test_that("ADV-D2: unknown priors_by_kind entry aborts with gdpar_input_error", {
  expect_error(
    gdpar_prior(priors_by_kind = list(not_a_kind = "normal(0, 1)")),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_prior(priors_by_kind = list(mu = 1)),  # not a character scalar
    class = "gdpar_input_error"
  )
})

test_that("ADV-D3: content-degenerate prior syntax passes R-side (Stan-side reject) -- documented boundary", {
  # gdpar_prior validates *structure* (non-empty character scalar, known
  # canonical kind), not the numerical *content* of the Stan syntax.
  # A negative scale is therefore accepted here and only rejected by the
  # Stan parser/sampler. This is the deliberate R-side / Stan-side
  # validation boundary; documented in block9_internal.md.
  pr <- gdpar_prior(sigma_y = "normal(0, -1)")
  expect_s3_class(pr, "gdpar_prior")
  expect_identical(pr$sigma_y, "normal(0, -1)")
})

# ---- Section E: HOM (homogeneity) latent stratification -----------------

test_that("ADV-E1: latent two-group stratified data does not abort the pre-fit checks (HOM is a fit-time assumption)", {
  # Construct data with a latent grouping that violates cross-group
  # homogeneity of the slope. The HOM assumption is a distributional /
  # fit-time regularity condition, NOT a pre-fit algebraic invariant, so
  # gdpar_check_identifiability is expected to PASS (the design is full
  # rank) -- the violation can only surface post-fit (e.g. via the
  # contraction diagnostic). This test documents that boundary.
  set.seed(123)
  n <- 300
  g <- rep(c(0L, 1L), each = n / 2L)
  x1 <- rnorm(n)
  # opposite-sign slopes per latent group: a single pooled slope is
  # mis-specified, but the design matrix is still well conditioned.
  y <- ifelse(g == 0L, 2 * x1, -2 * x1) + rnorm(n, sd = 0.5)
  df <- data.frame(x1 = x1, y = y)
  spec <- amm_spec(a = ~ x1)
  rep <- gdpar_check_identifiability(spec, df, formula_rhs = ~ x1)
  expect_true(rep$passed)
  expect_lt(rep$condition_number, 1e6)
})

.make_mv_data <- function(n = 20L, seed = 1L) {
  set.seed(seed)
  data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
}

.materialize_spec_W <- function(spec) {
  if (!is.null(spec$W)) {
    spec$W <- materialize_W_basis(spec$W, p = spec$p)
  }
  spec
}

test_that("assemble_stan_data() multi: gaussian p=2 produces expected shape and field set", {
  d <- .make_mv_data()
  spec <- .materialize_spec_W(amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = ~ x1),
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1", "x2")
  ))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$n, 20L)
  expect_equal(sd$p, 2L)
  expect_equal(sd$family_id, 1L)
  expect_equal(sd$use_a, 1L)
  expect_equal(sd$use_b, 1L)
  expect_equal(sd$use_W, 1L)
  expect_equal(sd$J_a_per_k, c(2L, 2L))
  expect_equal(sd$J_b_per_k, c(1L, 1L))
  expect_equal(sd$J_a_max, 2L)
  expect_equal(sd$J_b_max, 1L)
  expect_length(sd$Z_a, 2L)
  expect_equal(dim(sd$Z_a[[1L]]), c(20L, 2L))
  expect_equal(dim(sd$Z_b[[1L]]), c(20L, 1L))
  expect_equal(dim(sd$y_real), c(20L, 2L))
  expect_equal(dim(sd$y_int), c(20L, 2L))
  expect_equal(sd$dim_W, 4L)
  expect_equal(sd$W_per_k_dim, 2L)
  expect_equal(sd$d, 2L)
  expect_equal(sd$theta_anchor, c(0, 0))
  expect_equal(sd$use_dispersion_y, 1L)
  expect_equal(sd$use_dispersion_phi, 0L)
  expect_equal(sd$W_type_id, 1L)
  expect_equal(sd$W_n_knots_full, 0L)
  expect_equal(length(sd$W_knots_full), 0L)
  expect_equal(sd$W_degree, 2L)
})

test_that("assemble_stan_data() multi: per-k overrides yield padded Z_a with zeros in tail", {
  d <- .make_mv_data()
  dims <- override(
    dimwise(a = ~ x1, b = NULL),
    k = 2L, a = ~ x1 + x2 + x3, b = NULL
  )
  spec <- amm_spec(p = 2L, dims = dims)
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2 + x3)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$J_a_per_k, c(1L, 3L))
  expect_equal(sd$J_a_max, 3L)
  expect_equal(dim(sd$Z_a[[1L]]), c(20L, 3L))
  expect_equal(dim(sd$Z_a[[2L]]), c(20L, 3L))
  expect_equal(unname(sd$Z_a[[1L]][, 2:3]),
               matrix(0, nrow = 20L, ncol = 2L))
  expect_equal(unname(sd$Z_a[[1L]][, 1L]),
               unname(design$Z_a_list[[1L]][, 1L]))
})

test_that("assemble_stan_data() multi: poisson p=3 routes outcomes to y_int and zeros y_real", {
  d <- .make_mv_data(n = 15L)
  spec <- amm_spec(p = 3L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("poisson", p = 3L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(sample.int(20L, 15L * 3L, replace = TRUE) - 1L,
              nrow = 15L, ncol = 3L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$family_id, 2L)
  expect_equal(sd$use_dispersion_y, 0L)
  expect_equal(sd$use_dispersion_phi, 0L)
  expect_equal(storage.mode(sd$y_int), "integer")
  expect_equal(sd$y_int, y)
  expect_equal(sd$y_real, matrix(0, nrow = 15L, ncol = 3L))
})

test_that("assemble_stan_data() multi: neg_binomial_2 uses dispersion_phi flag", {
  d <- .make_mv_data(n = 15L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("neg_binomial_2", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(sample.int(5L, 15L * 2L, replace = TRUE) - 1L,
              nrow = 15L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$family_id, 3L)
  expect_equal(sd$use_dispersion_y, 0L)
  expect_equal(sd$use_dispersion_phi, 1L)
})

test_that("assemble_stan_data() multi: bernoulli routes to y_int with binary content", {
  d <- .make_mv_data(n = 12L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("bernoulli", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(sample.int(2L, 12L * 2L, replace = TRUE) - 1L,
              nrow = 12L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$family_id, 4L)
  expect_equal(sd$use_dispersion_y, 0L)
  expect_equal(sd$use_dispersion_phi, 0L)
  expect_equal(sd$y_int, y)
})

test_that("assemble_stan_data() multi: no W -> dim_W and d are zero", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$use_W, 0L)
  expect_equal(sd$dim_W, 0L)
  expect_equal(sd$W_per_k_dim, 0L)
  expect_equal(sd$d, 0L)
  expect_equal(dim(sd$X), c(20L, 0L))
})

test_that("assemble_stan_data() multi: theta_anchor as vector of length p is preserved", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = c(1, -1))
  expect_equal(sd$theta_anchor, c(1, -1))
})

test_that("assemble_stan_data() multi: theta_anchor scalar is broadcast to length p", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 3L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 3L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 3L), nrow = 20L, ncol = 3L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0.5)
  expect_equal(sd$theta_anchor, c(0.5, 0.5, 0.5))
})

test_that("assemble_stan_data() multi: theta_anchor wrong length aborts", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  expect_error(
    assemble_stan_data(design, fam, spec, y, theta_anchor = c(0, 0, 0)),
    class = "gdpar_input_error"
  )
})

test_that("assemble_stan_data() multi: univariate family with amm$p > 1 aborts", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam_uni <- gdpar_family("gaussian")
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  expect_error(
    assemble_stan_data(design, fam_uni, spec, y, theta_anchor = 0),
    class = "gdpar_input_error"
  )
})

test_that("assemble_stan_data() multi: family p mismatch with amm p aborts", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 3L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  expect_error(
    assemble_stan_data(design, fam, spec, y, theta_anchor = 0),
    class = "gdpar_input_error"
  )
})

test_that("assemble_stan_data() multi: y as vector instead of matrix aborts", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  expect_error(
    assemble_stan_data(design, fam, spec, rnorm(20L), theta_anchor = 0),
    class = "gdpar_input_error"
  )
})

test_that("assemble_stan_data() multi: y with wrong column count aborts", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y_bad <- matrix(rnorm(20L * 3L), nrow = 20L, ncol = 3L)
  expect_error(
    assemble_stan_data(design, fam, spec, y_bad, theta_anchor = 0),
    class = "gdpar_input_error"
  )
})

test_that("assemble_stan_data() multi: bspline W without boundary_knots aborts", {
  # Sub-phase 8.3.8 (2026-05-22): the legacy guard that rejected
  # bspline outright was removed; bspline is now supported via Cox-de
  # Boor recursion (D2 = I) provided the user supplies boundary_knots.
  # Without boundary_knots, .gdpar_resolve_W_stan_data aborts with
  # gdpar_input_error pointing the user at the canonical fix.
  d <- .make_mv_data()
  spec <- .materialize_spec_W(amm_spec(
    p = 2L, dims = dimwise(a = ~ x1, b = NULL),
    W = W_basis(type = "bspline", degree = 3, df = 4),
    x_vars = "x1"
  ))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  expect_error(
    assemble_stan_data(design, fam, spec, y, theta_anchor = 0),
    class = "gdpar_input_error"
  )
})

test_that("assemble_stan_data() multi: bspline W with boundary_knots succeeds", {
  # Sub-phase 8.3.8 (2026-05-22): bspline + p > 1 now lives end-to-end
  # in amm_distrib_multi.stan via apply_W_basis_diff. The resolver
  # validates that the projected theta_ref range fits inside
  # boundary_knots (D4 = delta+alpha); pass wide boundary_knots here.
  d <- .make_mv_data()
  spec <- .materialize_spec_W(amm_spec(
    p = 2L, dims = dimwise(a = ~ x1, b = NULL),
    W = W_basis(type = "bspline", degree = 3, df = 4,
                boundary_knots = c(-20, 20)),
    x_vars = "x1"
  ))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$W_type_id, 2L)
  expect_equal(sd$W_degree, 3L)
  # df = 4, degree = 3 -> n_int = 1; augmented knot vector length =
  # n_int + 2*(degree+1) = 1 + 8 = 9.
  expect_equal(sd$W_n_knots_full, 9L)
  expect_equal(sd$dim_W, 8L)
  expect_equal(sd$W_per_k_dim, 4L)
})

test_that("assemble_stan_data() multi: scalar path (p=1L) untouched by dispatcher", {
  d <- .make_mv_data()
  spec <- amm_spec(a = ~ x1, b = NULL)
  fam <- gdpar_family("gaussian")
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- rnorm(20L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_named(sd, c("n", "family_id", "use_a", "use_b", "use_W",
                      "J_a", "J_b", "dim_W", "d", "Z_a", "Z_b", "X",
                      "y_real", "y_int", "theta_anchor",
                      "use_dispersion_y", "use_dispersion_phi",
                      "W_type_id", "W_n_knots_full", "W_knots_full",
                      "W_degree",
                      "use_groups", "J_groups", "group_id"))
  expect_false("p" %in% names(sd))
  expect_false("J_a_per_k" %in% names(sd))
  expect_equal(sd$use_groups, 0L)
  expect_equal(sd$J_groups, 1L)
  expect_equal(length(sd$group_id), 20L)
  expect_true(all(sd$group_id == 1L))
})

test_that("generate_stan_code_multi() with default prior substitutes all placeholders", {
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_match(src, "data\\s*\\{")
  expect_match(src, "parameters\\s*\\{")
  expect_match(src, "transformed\\s+parameters\\s*\\{")
  expect_match(src, "model\\s*\\{")
  expect_match(src, "generated\\s+quantities\\s*\\{")
})

test_that("generate_stan_code_multi() NCP uses sigma_a[k] (per-coord scale), not sigma_a[1]", {
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior, cp_a = FALSE)
  expect_match(src, "* sigma_a[k]", fixed = TRUE)
  expect_false(grepl(" * sigma_a[1]", src, fixed = TRUE))
})

test_that("generate_stan_code_multi() NCP W uses sigma_W[1] (W is globally shared)", {
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior, cp_W = FALSE)
  expect_match(src, "* sigma_W[1]", fixed = TRUE)
})

test_that("generate_stan_code_multi() CP W swaps prior to normal(0, sigma_W[1]) and drops scale at use site", {
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior, cp_W = TRUE)
  expect_match(src, "normal(0, sigma_W[1])", fixed = TRUE)
})

test_that("generate_stan_code_multi() accepts cp_a = TRUE (Phase H.2 unblock)", {
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior, cp_a = TRUE)
  expect_type(src, "character")
  expect_match(src, "segment\\(a_raw")
})

test_that("generate_stan_code_multi() inserts custom prior strings verbatim", {
  prior <- gdpar_prior(theta_ref = "normal(1.5, 3)",
                       sigma_y   = "exponential(1)")
  src <- generate_stan_code_multi(prior)
  expect_match(src, "normal(1.5, 3)", fixed = TRUE)
  expect_match(src, "exponential(1)", fixed = TRUE)
})

test_that("generate_stan_code_multi() result is syntactically valid Stan", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior)
  tmpfile <- tempfile(fileext = ".stan")
  writeLines(src, tmpfile)
  mod <- cmdstanr::cmdstan_model(tmpfile, compile = FALSE)
  expect_true(mod$check_syntax(pedantic = FALSE, quiet = TRUE))
})

test_that("assemble_stan_data() multi: all-zero J_a_per_k yields J_a_max == 0 and empty Z_a matrices", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = NULL, b = ~ x1))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  y <- matrix(rnorm(20L * 2L), nrow = 20L, ncol = 2L)
  sd <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(sd$use_a, 0L)
  expect_equal(sd$J_a_max, 0L)
  expect_equal(sd$J_a_per_k, c(0L, 0L))
  for (k in seq_len(2L)) {
    expect_equal(dim(sd$Z_a[[k]]), c(20L, 0L))
  }
})

test_that("legacy K=1 FB templates relocated to _canonical_pieces/ under Sub-sub-fase 9.3.a", {
  # Sub-sub-fase 9.3.a (Bloque 9, Sesion B9.3, 2026-05-27) canonized B.iv
  # lateral: the legacy K=1 FB templates amm_main.stan (p=1) and
  # amm_distrib_multi.stan (p>=1) are relocated to
  # inst/stan/_canonical_pieces/ as amm_canonical_p1.stan and
  # amm_canonical_pmulti.stan respectively (byte-identical copies); the
  # legacy paths at the top-level of inst/stan/ are removed. The
  # canonical dispatcher .gdpar_emit_canonical_stan() reads from the
  # _canonical_pieces/ location. Sub-fase 8.2(b) preservation of the
  # amm_main_multi.stan removal is maintained.
  for (piece in c("amm_canonical_p1.stan", "amm_canonical_pmulti.stan")) {
    installed <- system.file("stan", "_canonical_pieces", piece,
                             package = "gdpar")
    source_path <- file.path("..", "..", "inst", "stan",
                             "_canonical_pieces", piece)
    expect_true(nzchar(installed) || file.exists(source_path),
                info = sprintf("canonical piece '%s' must exist", piece))
  }
  for (legacy in c("amm_main.stan", "amm_distrib_multi.stan",
                   "amm_main_multi.stan")) {
    old_installed <- system.file("stan", legacy, package = "gdpar")
    old_source <- file.path("..", "..", "inst", "stan", legacy)
    expect_false(nzchar(old_installed) && file.exists(old_installed),
                 info = sprintf("legacy template '%s' must be removed", legacy))
    expect_false(file.exists(old_source),
                 info = sprintf("legacy source '%s' must be removed", legacy))
  }
})

test_that(".gdpar_n_params_individual returns 1 for every built-in family (K=1 default in B8.0)", {
  for (name in c("gaussian", "poisson", "neg_binomial_2", "bernoulli")) {
    fam <- gdpar_family(name)
    expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  }
  fam_mv <- gdpar_family_multi("gaussian", p = 2L)
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_mv), 1L)
})

test_that(".gdpar_guard_multiparam_multivariate aborts when K>1 and p>1, passes otherwise", {
  fake_K2 <- gdpar_family("gaussian")
  fake_K2$param_specs[[2L]]$scope <- "per_group"
  expect_equal(gdpar:::.gdpar_n_params_individual(fake_K2), 2L)
  expect_error(
    gdpar:::.gdpar_guard_multiparam_multivariate(fake_K2, p = 2L),
    class = "gdpar_unsupported_feature_error"
  )
  expect_silent(
    gdpar:::.gdpar_guard_multiparam_multivariate(fake_K2, p = 1L)
  )
  expect_silent(
    gdpar:::.gdpar_guard_multiparam_multivariate(gdpar_family("gaussian"),
                                                 p = 2L)
  )
})

# Sub-phase 8.3.3 Unit 3 (2026-05-20): structural tests for the
# K-individual path (.build_amm_design_K, .assemble_stan_data_K,
# generate_a_blocks_K, generate_stan_code_K, resolve_anchor_K,
# gdpar() dispatch on K > 1 with p = 1). Stan compilation paths are
# gated by env vars and live in companion files.

.make_K_data <- function(n = 40L, seed = 7L) {
  set.seed(seed)
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    z = rnorm(n),
    w = rnorm(n)
  )
}

# .build_amm_design_K -------------------------------------------------------

test_that(".build_amm_design_K returns K, slot_names, per-slot lists", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x),
    sigma = amm_spec(a = ~ z)
  )
  out <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  expect_named(out, c("K", "slot_names", "Z_a_k_list", "Z_b_k_list",
                       "X", "Z_a_k_means_list", "Z_b_k_means_list",
                       "X_means", "X_sds",
                       "Z_a_k_names_list", "Z_b_k_names_list",
                       "X_names"))
  expect_equal(out$K, 2L)
  expect_equal(out$slot_names, c("mu", "sigma"))
  expect_length(out$Z_a_k_list, 2L)
  expect_length(out$Z_b_k_list, 2L)
  expect_equal(ncol(out$Z_a_k_list[[1L]]), 1L)
  expect_equal(ncol(out$Z_a_k_list[[2L]]), 1L)
  expect_equal(ncol(out$Z_b_k_list[[1L]]), 0L)
  expect_equal(ncol(out$Z_b_k_list[[2L]]), 0L)
  expect_equal(nrow(out$Z_a_k_list[[1L]]), 40L)
})

test_that(".build_amm_design_K supports ragged column counts across slots", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x + z, b = ~ w),
    sigma = amm_spec(a = ~ x)
  )
  out <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z + w)
  expect_equal(ncol(out$Z_a_k_list[[1L]]), 2L)
  expect_equal(ncol(out$Z_a_k_list[[2L]]), 1L)
  expect_equal(ncol(out$Z_b_k_list[[1L]]), 1L)
  expect_equal(ncol(out$Z_b_k_list[[2L]]), 0L)
})

test_that(".build_amm_design_K column-centers Z_a_k and Z_b_k per slot", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x + z),
    sigma = amm_spec(a = ~ z, b = ~ w)
  )
  out <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z + w)
  for (k in seq_len(out$K)) {
    if (ncol(out$Z_a_k_list[[k]]) > 0L) {
      cm <- colMeans(out$Z_a_k_list[[k]])
      expect_equal(unname(cm), rep(0, length(cm)), tolerance = 1e-12)
    }
    if (ncol(out$Z_b_k_list[[k]]) > 0L) {
      cm <- colMeans(out$Z_b_k_list[[k]])
      expect_equal(unname(cm), rep(0, length(cm)), tolerance = 1e-12)
    }
  }
})

test_that(".build_amm_design_K shares a single X across slots when W is active", {
  d <- .make_K_data()
  W_obj <- W_basis(type = "polynomial", degree = 2L, p = 1L)
  amm_list <- list(
    mu    = amm_spec(a = ~ x, W = W_obj),
    sigma = amm_spec(a = ~ z, W = W_obj)
  )
  out <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  expect_true(is.matrix(out$X))
  expect_equal(nrow(out$X), 40L)
  expect_true(ncol(out$X) >= 2L)
  expect_true(setequal(out$X_names, c("x", "z")))
})

test_that(".build_amm_design_K returns empty X when no slot declares W", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x),
    sigma = amm_spec(a = ~ z)
  )
  out <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  expect_equal(ncol(out$X), 0L)
  expect_equal(length(out$X_names), 0L)
})

test_that(".build_amm_design_K aborts when any slot carries p > 1", {
  d <- .make_K_data()
  amm_p2 <- amm_spec(p = 2L, dims = dimwise(a = ~ x, b = NULL))
  amm_list <- list(
    mu    = amm_p2,
    sigma = amm_spec(a = ~ z)
  )
  expect_error(
    .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".build_amm_design_K aborts when K = 1", {
  d <- .make_K_data()
  expect_error(
    .build_amm_design_K(list(mu = amm_spec(a = ~ x)), d,
                         formula_rhs = ~ x),
    class = "gdpar_internal_error"
  )
})

test_that(".build_amm_design_K aborts on missing values in needed covariates", {
  d <- .make_K_data()
  d$x[3L] <- NA
  amm_list <- list(
    mu    = amm_spec(a = ~ x),
    sigma = amm_spec(a = ~ z)
  )
  expect_error(
    .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z),
    class = "gdpar_input_error"
  )
})

test_that(".build_amm_design_K aborts on constant covariate in W block", {
  d <- .make_K_data()
  d$x <- 0  # constant
  W_obj <- W_basis(type = "polynomial", degree = 2L, p = 1L)
  amm_list <- list(
    mu    = amm_spec(a = ~ x, W = W_obj),
    sigma = amm_spec(a = ~ z, W = W_obj)
  )
  expect_error(
    .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z),
    class = "gdpar_input_error"
  )
})

# .assemble_stan_data_K ----------------------------------------------------

test_that(".assemble_stan_data_K returns the expected Stan keys for K=2 Gaussian", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x),
    sigma = amm_spec(a = ~ z)
  )
  design_K <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  fam <- .gdpar_promote_scope_per_observation(gdpar_family("gaussian"),
                                              c("mu", "sigma"))
  sd <- .assemble_stan_data_K(design_K, fam, amm_list, d$y,
                               c(mu = 0, sigma = 0))
  expect_named(sd,
               c("n", "K", "family_id_k", "inv_link_id_per_slot",
                 "use_a_k", "use_b_k", "use_W",
                 "J_a_per_k", "J_b_per_k", "J_a_max", "J_b_max",
                 "Z_a_k", "Z_b_k", "dim_W", "d", "X",
                 "y_real", "y_int", "theta_anchor_K",
                 "use_dispersion_y_k", "use_dispersion_phi_k",
                 "W_type_id", "W_n_knots_full", "W_knots_full",
                 "W_degree", "use_groups", "J_groups", "group_id"))
  expect_equal(sd$n, 40L)
  expect_equal(sd$K, 2L)
  expect_equal(sd$family_id_k, c(1L, 1L))  # gaussian, homogeneous
  expect_equal(sd$use_a_k, c(1L, 1L))
  expect_equal(sd$use_b_k, c(0L, 0L))
  expect_equal(sd$use_W, 0L)
  expect_equal(sd$J_a_per_k, c(1L, 1L))
  expect_equal(sd$J_a_max, 1L)
  expect_equal(length(sd$Z_a_k), 2L)
  expect_equal(dim(sd$Z_a_k[[1L]]), c(40L, 1L))
  expect_equal(sd$theta_anchor_K, c(0, 0))
  expect_equal(sd$use_dispersion_y_k, c(0L, 0L))
  expect_equal(sd$use_dispersion_phi_k, c(0L, 0L))
})

test_that(".assemble_stan_data_K pads ragged Z_*_k to J_*_max", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x + z),  # J_a = 2
    sigma = amm_spec(a = ~ x)        # J_a = 1
  )
  design_K <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  fam <- .gdpar_promote_scope_per_observation(gdpar_family("gaussian"),
                                              c("mu", "sigma"))
  sd <- .assemble_stan_data_K(design_K, fam, amm_list, d$y, c(0, 0))
  expect_equal(sd$J_a_per_k, c(2L, 1L))
  expect_equal(sd$J_a_max, 2L)
  expect_equal(ncol(sd$Z_a_k[[1L]]), 2L)
  expect_equal(ncol(sd$Z_a_k[[2L]]), 2L)  # padded
  # padded column on slot 2 is exactly zero
  expect_true(all(sd$Z_a_k[[2L]][, 2L] == 0))
})

test_that(".assemble_stan_data_K emits homogeneous family_id_k for neg_binomial_2 K=2", {
  d <- .make_K_data()
  d$y <- rpois(nrow(d), lambda = 3)
  amm_list <- list(
    mu  = amm_spec(a = ~ x),
    phi = amm_spec(a = ~ z)
  )
  design_K <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  fam <- .gdpar_promote_scope_per_observation(
    gdpar_family("neg_binomial_2"), c("mu", "phi")
  )
  sd <- .assemble_stan_data_K(design_K, fam, amm_list, d$y,
                               c(mu = 0, phi = 0))
  expect_equal(sd$family_id_k, c(3L, 3L))
  expect_true(all(sd$y_int == as.integer(d$y)))
  expect_true(all(sd$y_real == 0))
})

test_that(".assemble_stan_data_K aborts when family is gdpar_family_multi", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x),
    sigma = amm_spec(a = ~ z)
  )
  design_K <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  fam <- gdpar_family_multi(gdpar_family("gaussian"), p = 2L)
  expect_error(
    .assemble_stan_data_K(design_K, fam, amm_list, d$y, c(0, 0)),
    class = "gdpar_internal_error"
  )
})

test_that(".assemble_stan_data_K aborts on un-wired family (stan_id not in {1, 3})", {
  d <- .make_K_data()
  amm_list <- list(
    mu    = amm_spec(a = ~ x),
    sigma = amm_spec(a = ~ z)
  )
  design_K <- .build_amm_design_K(amm_list, d, formula_rhs = ~ x + z)
  # Forge a family with stan_id = 2 (poisson) just to test the guard
  fam <- gdpar_family("poisson")
  fam$param_specs <- list(
    list(name = "mu",    scope = "per_observation"),
    list(name = "sigma", scope = "per_observation")
  )
  expect_error(
    .assemble_stan_data_K(design_K, fam, amm_list, d$y, c(0, 0)),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".assemble_stan_data_K aborts on K = 1", {
  d <- .make_K_data()
  amm_list <- list(mu = amm_spec(a = ~ x))
  design_K_fake <- list(K = 1L, slot_names = "mu",
                        Z_a_k_list = list(matrix(0, 40, 1)),
                        Z_b_k_list = list(matrix(0, 40, 0)),
                        X = matrix(0, 40, 0))
  fam <- .gdpar_promote_scope_per_observation(gdpar_family("gaussian"),
                                              "mu")
  expect_error(
    .assemble_stan_data_K(design_K_fake, fam, amm_list, d$y, 0),
    class = "gdpar_internal_error"
  )
})

# generate_a_blocks_K ------------------------------------------------------

test_that("generate_a_blocks_K uniform-NCP path emits empty data_decl", {
  out <- generate_a_blocks_K(c(FALSE, FALSE))
  expect_named(out, c("data_decl", "tp_block", "model_block"))
  expect_equal(out$data_decl, "")
  expect_match(out$model_block, "a_raw ~ normal\\(0, 1\\)")
  # Option A (D96): sigma_a_k is compacted and indexed via sigma_a_idx[k].
  expect_match(out$tp_block, "sigma_a_k\\[sigma_a_idx\\[k\\]\\]")
})

test_that("generate_a_blocks_K uniform-CP path uses segment-based prior", {
  out <- generate_a_blocks_K(c(TRUE, TRUE))
  expect_equal(out$data_decl, "")
  expect_match(out$model_block, "segment\\(a_raw")
  expect_match(out$model_block, "normal\\(0, sigma_a_k\\[sigma_a_idx\\[k\\]\\]\\)")
  expect_false(grepl("a_raw\\[a_raw_offset", out$tp_block) &&
                grepl("sigma_a_k\\[sigma_a_idx\\[k\\]\\]", out$tp_block))
})

test_that("generate_a_blocks_K mixed path declares cp_a_per_K_data", {
  out <- generate_a_blocks_K(c(TRUE, FALSE))
  expect_match(out$data_decl, "cp_a_per_K_data")
  expect_match(out$tp_block, "cp_a_per_K_data\\[k\\] == 1")
  expect_match(out$model_block, "cp_a_per_K_data\\[k\\] == 1")
})

# generate_stan_code_K -----------------------------------------------------

test_that("generate_stan_code_K substitutes all placeholders for default prior", {
  src <- generate_stan_code_K(gdpar_prior())
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(nchar(src) > 500L)
  expect_match(src, "amm_distrib_K", fixed = FALSE)
})

test_that("generate_stan_code_K honors cp_a_per_K mixed pattern", {
  src <- generate_stan_code_K(gdpar_prior(), cp_a_per_K = c(TRUE, FALSE))
  expect_match(src, "cp_a_per_K_data")
})

test_that("generate_stan_code_K aborts on malformed cp_a_per_K argument", {
  expect_error(
    generate_stan_code_K(gdpar_prior(), cp_a_per_K = c(NA, TRUE)),
    class = "gdpar_input_error"
  )
})

# resolve_anchor_K ---------------------------------------------------------

test_that("resolve_anchor_K broadcasts a scalar to length K", {
  out <- resolve_anchor_K(0, gdpar_family("gaussian"), rnorm(20),
                          gdpar_prior(), c("mu", "sigma"),
                          verbose = FALSE)
  expect_equal(out, c(mu = 0, sigma = 0))
})

test_that("resolve_anchor_K accepts an unnamed numeric vector of length K", {
  out <- resolve_anchor_K(c(1.5, -2), gdpar_family("gaussian"),
                          rnorm(20), gdpar_prior(),
                          c("mu", "sigma"), verbose = FALSE)
  expect_equal(unname(out), c(1.5, -2))
  expect_equal(names(out), c("mu", "sigma"))
})

test_that("resolve_anchor_K reorders a named numeric vector by slot names", {
  out <- resolve_anchor_K(c(sigma = 3, mu = 7),
                          gdpar_family("gaussian"),
                          rnorm(20), gdpar_prior(),
                          c("mu", "sigma"), verbose = FALSE)
  expect_equal(out[["mu"]], 7)
  expect_equal(out[["sigma"]], 3)
})

test_that("resolve_anchor_K rejects a named numeric vector with non-matching names", {
  expect_error(
    resolve_anchor_K(c(alpha = 1, beta = 2),
                     gdpar_family("gaussian"), rnorm(20),
                     gdpar_prior(), c("mu", "sigma"),
                     verbose = FALSE),
    class = "gdpar_input_error"
  )
})

test_that("resolve_anchor_K maps 'prior_mean' to a vector of zeros", {
  out <- resolve_anchor_K("prior_mean", gdpar_family("gaussian"),
                          rnorm(20), gdpar_prior(),
                          c("mu", "sigma"), verbose = FALSE)
  expect_equal(unname(out), c(0, 0))
})

test_that("resolve_anchor_K 'empirical_y' places link(mean(y)) on the location slot", {
  set.seed(13L)
  y <- rnorm(50L, mean = 2.5)
  out <- resolve_anchor_K("empirical_y", gdpar_family("gaussian"),
                          y, gdpar_prior(),
                          c("mu", "sigma"), verbose = FALSE)
  expect_equal(unname(out[["mu"]]), mean(y))
  expect_equal(unname(out[["sigma"]]), 0)
})

test_that("resolve_anchor_K rejects nonsense inputs with gdpar_input_error", {
  expect_error(
    resolve_anchor_K(list(mu = 1, sigma = 2),
                     gdpar_family("gaussian"), rnorm(20),
                     gdpar_prior(), c("mu", "sigma"),
                     verbose = FALSE),
    class = "gdpar_input_error"
  )
})

# gdpar() dispatch on K > 1 ------------------------------------------------

test_that("gdpar() with K=2 formula set + Gaussian routes to .gdpar_K (no Stan)", {
  # We verify the dispatch path without compiling Stan by mocking
  # cmdstanr::cmdstan_model to throw a sentinel error and catching it.
  d <- .make_K_data()
  skip_if_not_installed("cmdstanr")
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  sentinel <- "gdpar_K_dispatch_sentinel"
  hit <- tryCatch(
    testthat::with_mocked_bindings(
      cmdstan_model = function(...) stop(sentinel, call. = FALSE),
      .package = "cmdstanr",
      code = {
        gdpar(fs, data = d, family = gdpar_family("gaussian"),
              skip_id_check = TRUE, verbose = FALSE, seed = 1L,
              chains = 1L, iter_warmup = 10L, iter_sampling = 10L)
      }
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(sentinel, hit, fixed = TRUE))
})

test_that("gdpar() with K=2 + non-eligible slot aborts via scope promotion", {
  d <- .make_K_data()
  fs <- gdpar_bf(y ~ a(x), nonexistent ~ a(z))
  expect_error(
    gdpar(fs, data = d, family = gdpar_family("gaussian"),
          skip_id_check = TRUE, verbose = FALSE),
    class = "gdpar_input_error"
  )
})

test_that("gdpar() with K=2 + gdpar_family_multi aborts (K>1+p>1 deferred)", {
  d <- .make_K_data()
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fam_multi <- gdpar_family_multi(gdpar_family("gaussian"), p = 2L)
  # The dispatch invokes .assemble_stan_data_K which aborts on
  # gdpar_family_multi inputs (internal_error per Unit 3 design).
  expect_error(
    gdpar(fs, data = d, family = fam_multi,
          skip_id_check = TRUE, verbose = FALSE),
    class = "gdpar_internal_error"
  )
})

test_that(".gdpar_promote_scope_per_observation promotes both K=2 slots", {
  fam <- gdpar_family("gaussian")
  fam2 <- .gdpar_promote_scope_per_observation(fam, c("mu", "sigma"))
  scopes <- vapply(fam2$param_specs, function(s) s$scope, character(1L))
  names(scopes) <- vapply(fam2$param_specs, function(s) s$name,
                          character(1L))
  expect_equal(scopes[["mu"]], "per_observation")
  expect_equal(scopes[["sigma"]], "per_observation")
})

# Backward-compat sanity: K=1 path still routes legacy ---------------------

test_that("gdpar() with K=1 formula set still routes to legacy K=1 path", {
  d <- .make_K_data()
  skip_if_not_installed("cmdstanr")
  fs <- gdpar_bf(y ~ a(x))
  sentinel <- "gdpar_K1_legacy_sentinel"
  hit <- tryCatch(
    testthat::with_mocked_bindings(
      cmdstan_model = function(...) stop(sentinel, call. = FALSE),
      .package = "cmdstanr",
      code = {
        gdpar(fs, data = d, family = gdpar_family("gaussian"),
              skip_id_check = TRUE, verbose = FALSE, seed = 1L,
              chains = 1L, iter_warmup = 10L, iter_sampling = 10L)
      }
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(sentinel, hit, fixed = TRUE))
})

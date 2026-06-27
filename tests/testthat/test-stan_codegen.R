test_that("generate_stan_code substitutes all placeholders", {
  pr <- gdpar_prior()
  src <- gdpar:::generate_stan_code(pr)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("normal(0, 2.5)", src, fixed = TRUE))
})

test_that("generate_stan_code default emits NCP form for a", {
  pr <- gdpar_prior()
  src <- gdpar:::generate_stan_code(pr)
  expect_true(grepl("a_coef[j] = a_raw[j] * sigma_a[1];",
                    src, fixed = TRUE))
  expect_true(grepl("a_raw ~ normal(0, 1);", src, fixed = TRUE))
  expect_false(grepl("a_coef[j] = a_raw[j];", src, fixed = TRUE))
  expect_false(grepl("a_raw ~ normal(0, sigma_a[1]);",
                     src, fixed = TRUE))
  # (C2) is enforced by column-wise centering of Z_a in R, not by a
  # hardcoded sum-to-zero on a_coef. Verify the redundant constraint
  # is NOT emitted in either parametrization.
  expect_false(grepl("a_coef[J_a] = -sum(a_raw) * sigma_a[1];",
                     src, fixed = TRUE))
  expect_false(grepl("a_coef[J_a] = -sum(a_raw);",
                     src, fixed = TRUE))
})

test_that("generate_stan_code with cp_a = TRUE switches to CP for a", {
  pr <- gdpar_prior()
  src <- gdpar:::generate_stan_code(pr, cp_a = TRUE)
  expect_true(grepl("a_coef[j] = a_raw[j];", src, fixed = TRUE))
  expect_true(grepl("a_raw ~ normal(0, sigma_a[1]);",
                    src, fixed = TRUE))
  expect_false(grepl("a_coef[j] = a_raw[j] * sigma_a[1];",
                     src, fixed = TRUE))
  expect_false(grepl("a_raw ~ normal(0, 1);", src, fixed = TRUE))
  expect_false(grepl("{{", src, fixed = TRUE))
  # Same redundant sum-to-zero must remain absent under CP.
  expect_false(grepl("a_coef[J_a] = -sum(a_raw);",
                     src, fixed = TRUE))
  expect_false(grepl("a_coef[J_a] = -sum(a_raw) * sigma_a[1];",
                     src, fixed = TRUE))
})

test_that("generate_stan_code default emits NCP form for W", {
  pr <- gdpar_prior()
  src <- gdpar:::generate_stan_code(pr)
  expect_true(grepl(
    "W_diff_x[jj] += basis_diff[k] * W_raw[k, jj] * sigma_W[1];",
    src, fixed = TRUE
  ))
  expect_true(grepl("to_vector(W_raw) ~ normal(0, 1);",
                    src, fixed = TRUE))
  expect_false(grepl("W_diff_x[jj] += basis_diff[k] * W_raw[k, jj];",
                     src, fixed = TRUE))
  expect_false(grepl("to_vector(W_raw) ~ normal(0, sigma_W[1]);",
                     src, fixed = TRUE))
})

test_that("generate_stan_code with cp_W = TRUE switches to CP for W", {
  pr <- gdpar_prior()
  src <- gdpar:::generate_stan_code(pr, cp_W = TRUE)
  expect_true(grepl("W_diff_x[jj] += basis_diff[k] * W_raw[k, jj];",
                    src, fixed = TRUE))
  expect_true(grepl("to_vector(W_raw) ~ normal(0, sigma_W[1]);",
                    src, fixed = TRUE))
  expect_false(grepl(
    "W_diff_x[jj] += basis_diff[k] * W_raw[k, jj] * sigma_W[1];",
    src, fixed = TRUE
  ))
  expect_false(grepl("to_vector(W_raw) ~ normal(0, 1);",
                     src, fixed = TRUE))
  expect_false(grepl("{{", src, fixed = TRUE))
})

test_that("generate_stan_code composes cp_a and cp_W independently", {
  pr <- gdpar_prior()
  src <- gdpar:::generate_stan_code(pr, cp_a = TRUE, cp_W = TRUE)
  expect_true(grepl("a_coef[j] = a_raw[j];", src, fixed = TRUE))
  expect_true(grepl("a_raw ~ normal(0, sigma_a[1]);",
                    src, fixed = TRUE))
  expect_true(grepl("W_diff_x[jj] += basis_diff[k] * W_raw[k, jj];",
                    src, fixed = TRUE))
  expect_true(grepl("to_vector(W_raw) ~ normal(0, sigma_W[1]);",
                    src, fixed = TRUE))
  expect_false(grepl("a_coef[j] = a_raw[j] * sigma_a[1];",
                     src, fixed = TRUE))
  expect_false(grepl(
    "W_diff_x[jj] += basis_diff[k] * W_raw[k, jj] * sigma_W[1];",
    src, fixed = TRUE
  ))
  expect_false(grepl("{{", src, fixed = TRUE))
})

test_that("generate_stan_code uses user-supplied prior strings", {
  pr <- gdpar_prior(theta_ref = "normal(1, 5)")
  src <- gdpar:::generate_stan_code(pr)
  expect_true(grepl("normal(1, 5)", src, fixed = TRUE))
})

test_that("write_stan_to_tempfile writes to tempdir and is reproducible", {
  src <- "// test stan code"
  p1 <- gdpar:::write_stan_to_tempfile(src)
  p2 <- gdpar:::write_stan_to_tempfile(src)
  expect_equal(p1, p2)
  expect_true(startsWith(p1, tempdir()))
})

test_that("assemble_stan_data builds expected structure", {
  df <- make_gaussian_data(n = 50, seed = 42)
  spec <- amm_spec(a = ~ x1 + x2)
  design <- gdpar:::build_amm_design(spec, df,
                                      formula_rhs = ~ x1 + x2)
  fam <- gdpar_family("gaussian")
  data_list <- gdpar:::assemble_stan_data(
    design, fam, spec, df$y, theta_anchor = 0
  )
  expect_equal(data_list$n, 50L)
  expect_equal(data_list$family_id, 1L)
  expect_equal(data_list$use_a, 1L)
  expect_equal(data_list$use_b, 0L)
  expect_equal(data_list$use_W, 0L)
  expect_equal(data_list$J_a, 2L)
  expect_equal(length(data_list$y_real), 50L)
})

test_that("assemble_stan_data rejects non-polynomial W in v0", {
  df <- make_gaussian_data(n = 50, seed = 42)
  spec <- amm_spec(
    W = W_basis(type = "user",
                basis_fn = function(theta) c(theta, theta^2),
                dim = 2),
    x_vars = "x1"
  )
  spec$W <- gdpar:::materialize_W_basis(spec$W, p = 1L)
  design <- gdpar:::build_amm_design(spec, df, formula_rhs = ~ x1 + x2)
  fam <- gdpar_family("gaussian")
  expect_error(
    gdpar:::assemble_stan_data(design, fam, spec, df$y,
                                theta_anchor = 0),
    class = "gdpar_unsupported_feature_error"
  )
})

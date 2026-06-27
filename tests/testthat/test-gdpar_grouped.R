make_grouped_gaussian_data <- function(n_per_group = 60L,
                                        J_groups = 3L,
                                        seed = 11L) {
  set.seed(seed)
  n <- n_per_group * J_groups
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  group <- factor(rep(letters[seq_len(J_groups)], each = n_per_group))
  theta_ref_g <- rnorm(J_groups, mean = 0, sd = 1.5)
  eta <- theta_ref_g[as.integer(group)] + 0.5 * x1 - 0.3 * x2
  y <- eta + rnorm(n, sd = 0.4)
  data.frame(x1 = x1, x2 = x2, group = group, y = y)
}

# ---- .resolve_group_argument: validation -----------------------------

test_that(".resolve_group_argument returns NULL when group is NULL", {
  d <- make_grouped_gaussian_data()
  expect_null(.resolve_group_argument(NULL, d, n = nrow(d),
                                       verbose = FALSE))
})

test_that(".resolve_group_argument rejects non-formula input", {
  d <- make_grouped_gaussian_data()
  expect_error(.resolve_group_argument("group", d, n = nrow(d),
                                        verbose = FALSE),
               class = "gdpar_input_error")
  expect_error(.resolve_group_argument(d$group, d, n = nrow(d),
                                        verbose = FALSE),
               class = "gdpar_input_error")
})

test_that(".resolve_group_argument rejects two-sided formula", {
  d <- make_grouped_gaussian_data()
  expect_error(.resolve_group_argument(y ~ group, d, n = nrow(d),
                                        verbose = FALSE),
               class = "gdpar_input_error")
})

test_that(".resolve_group_argument rejects formula with multiple vars", {
  d <- make_grouped_gaussian_data()
  d$group2 <- sample(letters[1:2], nrow(d), replace = TRUE)
  expect_error(.resolve_group_argument(~ group + group2, d,
                                        n = nrow(d), verbose = FALSE),
               class = "gdpar_input_error")
})

test_that(".resolve_group_argument rejects missing variable", {
  d <- make_grouped_gaussian_data()
  expect_error(.resolve_group_argument(~ does_not_exist, d,
                                        n = nrow(d), verbose = FALSE),
               class = "gdpar_input_error")
})

test_that(".resolve_group_argument rejects NA in grouping variable", {
  d <- make_grouped_gaussian_data()
  d$group[1] <- NA
  expect_error(.resolve_group_argument(~ group, d, n = nrow(d),
                                        verbose = FALSE),
               class = "gdpar_input_error")
})

test_that(".resolve_group_argument resolves to integer codes and levels", {
  d <- make_grouped_gaussian_data(n_per_group = 5L, J_groups = 3L)
  res <- .resolve_group_argument(~ group, d, n = nrow(d),
                                  verbose = FALSE)
  expect_type(res, "list")
  expect_equal(res$var_name, "group")
  expect_equal(res$levels, c("a", "b", "c"))
  expect_equal(length(res$group_id), nrow(d))
  expect_equal(sort(unique(res$group_id)), 1:3)
})

test_that(".resolve_group_argument warns when min(n_per_group) < 5", {
  d <- make_grouped_gaussian_data(n_per_group = 3L, J_groups = 3L)
  expect_warning(
    .resolve_group_argument(~ group, d, n = nrow(d), verbose = TRUE),
    class = "gdpar_grouping_warning"
  )
})

test_that(".resolve_group_argument is silent when min(n_per_group) >= 5", {
  d <- make_grouped_gaussian_data(n_per_group = 10L, J_groups = 3L)
  expect_silent(.resolve_group_argument(~ group, d, n = nrow(d),
                                         verbose = TRUE))
})

# ---- (C7) anti-aliasing ----------------------------------------------

test_that(".check_group_aliasing_c7 passes when Z_a has no group aliasing", {
  d <- make_grouped_gaussian_data(n_per_group = 20L, J_groups = 3L)
  spec <- amm_spec(a = ~ x1 + x2)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  group_id <- as.integer(d$group)
  expect_silent(.check_group_aliasing_c7(design, group_id, "group"))
})

test_that(".check_group_aliasing_c7 rejects factor(group) in a", {
  d <- make_grouped_gaussian_data(n_per_group = 20L, J_groups = 3L)
  spec <- amm_spec(a = ~ x1 + group)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1 + group)
  group_id <- as.integer(d$group)
  expect_error(.check_group_aliasing_c7(design, group_id, "group"),
               class = "gdpar_input_error")
})

test_that(".check_group_aliasing_c7 detects constant-per-group columns in b", {
  d <- make_grouped_gaussian_data(n_per_group = 20L, J_groups = 3L)
  d$z_const <- as.numeric(d$group)
  spec <- amm_spec(a = ~ x1, b = ~ z_const)
  design <- build_amm_design(spec, d,
                              formula_rhs = ~ x1 + z_const)
  group_id <- as.integer(d$group)
  expect_error(.check_group_aliasing_c7(design, group_id, "group"),
               class = "gdpar_input_error")
})

# ---- gdpar() integration: argument plumbing --------------------------

test_that("gdpar() accepts group = NULL silently (backward-compat path)", {
  d <- make_grouped_gaussian_data(n_per_group = 8L)
  spec <- amm_spec(a = ~ x1 + x2)
  skip_if_not_installed("cmdstanr")
  expect_silent({
    grp_info <- .resolve_group_argument(NULL, d, n = nrow(d),
                                         verbose = FALSE)
  })
  expect_null(grp_info)
})

test_that("gdpar() resolves group = ~ name into group_info structure", {
  d <- make_grouped_gaussian_data(n_per_group = 10L, J_groups = 3L)
  grp <- .resolve_group_argument(~ group, d, n = nrow(d),
                                  verbose = FALSE)
  expect_equal(grp$levels, c("a", "b", "c"))
  expect_equal(grp$var_name, "group")
})

# ---- Stan data assembly with group_id ---------------------------------

test_that("assemble_stan_data fills group_id from external arg", {
  d <- make_grouped_gaussian_data(n_per_group = 10L, J_groups = 3L)
  spec <- amm_spec(a = ~ x1)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  fam <- gdpar_family("gaussian")
  group_id <- as.integer(d$group)
  sd <- assemble_stan_data(design, fam, spec, d$y,
                            theta_anchor = 0,
                            group_id = group_id)
  expect_equal(sd$use_groups, 1L)
  expect_equal(sd$J_groups, 3L)
  expect_equal(sd$group_id, group_id)
})

test_that("assemble_stan_data backward-compat: group_id = NULL", {
  d <- make_grouped_gaussian_data(n_per_group = 10L, J_groups = 1L)
  spec <- amm_spec(a = ~ x1)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1)
  fam <- gdpar_family("gaussian")
  sd <- assemble_stan_data(design, fam, spec, d$y, theta_anchor = 0)
  expect_equal(sd$use_groups, 0L)
  expect_equal(sd$J_groups, 1L)
  expect_equal(length(sd$group_id), nrow(d))
  expect_true(all(sd$group_id == 1L))
})

# ---- generate_stan_code with sigma_theta_ref placeholder -------------

test_that("generate_stan_code substitutes the sigma_theta_ref placeholder", {
  pr <- gdpar_prior(sigma_theta_ref = "student_t(4, 0, 2)")
  src <- generate_stan_code(pr)
  expect_false(grepl("{{PRIOR_SIGMA_THETA_REF}}", src, fixed = TRUE))
  expect_true(grepl("student_t(4, 0, 2)", src, fixed = TRUE))
})

test_that("generate_stan_code_multi substitutes sigma_theta_ref placeholder", {
  pr <- gdpar_prior(sigma_theta_ref = "student_t(4, 0, 2)")
  src <- generate_stan_code_multi(pr)
  expect_false(grepl("{{PRIOR_SIGMA_THETA_REF}}", src, fixed = TRUE))
  expect_true(grepl("student_t(4, 0, 2)", src, fixed = TRUE))
})

# ---- gdpar_coef constructor under grouping ---------------------------

test_that("build_coef_theta_ref_df_grouped: J_groups = 1 collapses to (k,...) schema", {
  arr <- array(rnorm(40L * 1L * 2L), dim = c(40L, 1L, 2L))
  df <- build_coef_theta_ref_df_grouped(arr, J_groups = 1L, p = 2L)
  expect_named(df, c("k", "mean", "q05", "q50", "q95"))
  expect_equal(df$k, 1:2)
})

test_that("build_coef_theta_ref_df_grouped: J_groups > 1, p = 1 emits matrix path", {
  arr <- matrix(rnorm(40L * 3L), nrow = 40L, ncol = 3L)
  df <- build_coef_theta_ref_df_grouped(arr, J_groups = 3L, p = 1L)
  expect_named(df, c("g", "k", "mean", "q05", "q50", "q95"))
  expect_equal(nrow(df), 3L)
  expect_equal(df$g, 1:3)
  expect_true(all(df$k == 1L))
})

test_that("build_coef_theta_ref_df_grouped: J_groups > 1, p > 1 emits (g,k) schema", {
  arr <- array(rnorm(40L * 3L * 2L), dim = c(40L, 3L, 2L))
  df <- build_coef_theta_ref_df_grouped(arr, J_groups = 3L, p = 2L)
  expect_named(df, c("g", "k", "mean", "q05", "q50", "q95"))
  expect_equal(nrow(df), 6L)
  expect_equal(sort(unique(df$g)), 1:3)
  expect_equal(sort(unique(df$k)), 1:2)
})

test_that("new_gdpar_coef accepts mu_theta_ref and sigma_theta_ref slots", {
  arr <- array(rnorm(40L * 2L * 1L), dim = c(40L, 2L, 1L))
  theta_df <- build_coef_theta_ref_df_grouped(arr, J_groups = 2L,
                                               p = 1L)
  mu_df <- data.frame(k = 1L, mean = 0.1, q05 = -0.1, q50 = 0.1,
                      q95 = 0.3)
  sigma_df <- data.frame(k = 1L, mean = 1, q05 = 0.5, q50 = 1,
                         q95 = 1.5)
  obj <- new_gdpar_coef(
    theta_ref       = theta_df,
    p               = 1L,
    mu_theta_ref    = mu_df,
    sigma_theta_ref = sigma_df,
    J_groups        = 2L,
    group_levels    = c("a", "b")
  )
  expect_s3_class(obj, "gdpar_coef")
  expect_equal(obj$J_groups, 2L)
  expect_equal(obj$group_levels, c("a", "b"))
  expect_s3_class(obj$mu_theta_ref, "data.frame")
  expect_s3_class(obj$sigma_theta_ref, "data.frame")
})

test_that("new_gdpar_coef rejects hyperparam slots with wrong nrow", {
  arr <- array(rnorm(40L * 2L * 1L), dim = c(40L, 2L, 1L))
  theta_df <- build_coef_theta_ref_df_grouped(arr, J_groups = 2L,
                                               p = 1L)
  bad <- data.frame(k = 1:2, mean = c(0.1, 0.2), q05 = c(0, 0),
                    q50 = c(0.1, 0.2), q95 = c(0.3, 0.4))
  expect_error(
    new_gdpar_coef(theta_ref    = theta_df,
                   p            = 1L,
                   mu_theta_ref = bad,
                   J_groups     = 2L),
    class = "gdpar_internal_error"
  )
})

test_that("print.gdpar_coef shows hyperparam blocks when present", {
  arr <- array(rnorm(40L * 2L * 1L), dim = c(40L, 2L, 1L))
  theta_df <- build_coef_theta_ref_df_grouped(arr, J_groups = 2L,
                                               p = 1L)
  mu_df <- data.frame(k = 1L, mean = 0.1, q05 = -0.1, q50 = 0.1,
                      q95 = 0.3)
  sigma_df <- data.frame(k = 1L, mean = 1, q05 = 0.5, q50 = 1,
                         q95 = 1.5)
  obj <- new_gdpar_coef(
    theta_ref       = theta_df,
    p               = 1L,
    mu_theta_ref    = mu_df,
    sigma_theta_ref = sigma_df,
    J_groups        = 2L
  )
  out <- capture.output(print(obj))
  expect_true(any(grepl("mu_theta_ref", out)))
  expect_true(any(grepl("sigma_theta_ref", out)))
  expect_true(any(grepl("J_groups", out)))
})

# ---- End-to-end smoke fit (env-gated) --------------------------------

test_that("gdpar() end-to-end with group: fit, coef, predict cross-group", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          "Set GDPAR_F6_FIT=1 to run the grouped end-to-end smoke (compiles a Stan model; takes minutes).")
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  d <- make_grouped_gaussian_data(n_per_group = 40L, J_groups = 3L,
                                   seed = 13L)
  spec <- amm_spec(a = ~ x1 + x2)
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family("gaussian"),
    amm           = spec,
    data          = d,
    group         = ~ group,
    chains        = 2L,
    iter_warmup   = 200L,
    iter_sampling = 200L,
    refresh       = 0L,
    seed          = 13L,
    verbose       = FALSE
  )
  expect_s3_class(fit, "gdpar_fit")
  expect_false(is.null(fit$group_info))
  expect_equal(fit$group_info$var_name, "group")
  expect_equal(length(fit$group_info$levels), 3L)
  cf <- coef(fit)
  expect_s3_class(cf$theta_ref, "data.frame")
  expect_true("g" %in% names(cf$theta_ref))
  expect_equal(sort(unique(cf$theta_ref$g)), 1:3)
  expect_s3_class(cf$mu_theta_ref, "data.frame")
  expect_s3_class(cf$sigma_theta_ref, "data.frame")
  pred_in <- predict(fit, newdata = d, summary = "draws")
  expect_true(is.matrix(pred_in))
  d_new <- d[1:10L, ]
  d_new$group <- factor("z", levels = c(levels(d$group), "z"))
  expect_warning(
    pred_unseen <- predict(fit, newdata = d_new, summary = "draws"),
    class = "gdpar_predict_unseen_group_warning"
  )
  expect_equal(dim(pred_unseen)[2L], nrow(d_new))
})

test_that("gdpar() with group rejects factor(group) in a (C7 or id_check)", {
  d <- make_grouped_gaussian_data(n_per_group = 20L, J_groups = 3L)
  spec <- amm_spec(a = ~ x1 + group)
  skip_if_not_installed("cmdstanr")
  err <- tryCatch(
    gdpar(
      formula       = y ~ x1 + group,
      family        = gdpar_family("gaussian"),
      amm           = spec,
      data          = d,
      group         = ~ group,
      chains        = 1L,
      iter_warmup   = 50L,
      iter_sampling = 50L,
      refresh       = 0L,
      seed          = 17L,
      verbose       = FALSE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_true(any(c("gdpar_input_error", "gdpar_identifiability_error")
                  %in% class(err)))
})

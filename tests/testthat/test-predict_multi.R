test_that("gdpar() multi end-to-end p=2 gaussian: predict, summary, print, coef", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          "Set GDPAR_F6_FIT=1 to run the multi end-to-end smoke (compiles a Stan model; takes minutes).")
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  n <- 300L
  p <- 2L
  set.seed(42L)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  data_df <- data.frame(x1 = x1, x2 = x2)
  theta_ref_true <- c(0.5, -0.5)
  beta_a_true <- list(c(0.5, -0.5), c(-0.4, 0.4))
  sigma_true <- 0.3
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1_c + beta_a_true[[1]][2] * x2_c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1_c + beta_a_true[[2]][2] * x2_c
  )
  data_df$y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))

  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family_multi("gaussian", p = 2L),
    amm           = spec,
    data          = data_df,
    chains        = 2L,
    iter_warmup   = 200L,
    iter_sampling = 200L,
    refresh       = 0L,
    seed          = 42L,
    verbose       = FALSE
  )

  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$p, 2L)
  expect_s3_class(fit$family, "gdpar_family_multi")
  expect_lt(fit$diagnostics$rhat_max, 1.10)

  pred_draws <- predict(fit, summary = "draws")
  expect_true(is.array(pred_draws))
  expect_equal(length(dim(pred_draws)), 3L)
  expect_equal(dim(pred_draws)[2L], n)
  expect_equal(dim(pred_draws)[3L], p)
  expect_equal(dimnames(pred_draws)[[3L]], c("dim_1", "dim_2"))
  S <- dim(pred_draws)[1L]
  expect_equal(S, 400L)

  pred_newdata <- predict(fit, newdata = data_df, summary = "draws")
  expect_equal(dim(pred_newdata), dim(pred_draws))
  diff <- max(abs(pred_newdata - pred_draws))
  expect_lt(diff, 1e-6)

  pred_mean_se <- predict(fit, summary = "mean_se")
  expect_type(pred_mean_se, "list")
  expect_length(pred_mean_se, p)
  expect_named(pred_mean_se, c("dim_1", "dim_2"))
  for (k in seq_len(p)) {
    expect_s3_class(pred_mean_se[[k]], "data.frame")
    expect_equal(nrow(pred_mean_se[[k]]), n)
    expect_named(pred_mean_se[[k]], c("mean", "se"))
  }
  pred_q <- predict(fit, summary = "quantiles")
  expect_length(pred_q, p)
  for (k in seq_len(p)) {
    expect_s3_class(pred_q[[k]], "data.frame")
    expect_named(pred_q[[k]], c("q05", "q50", "q95"))
    expect_equal(nrow(pred_q[[k]]), n)
  }

  pred_resp <- predict(fit, type = "response", summary = "draws")
  expect_equal(pred_resp, pred_draws)

  summ <- summary(fit)
  expect_s3_class(summ, "data.frame")
  expect_true(all(c("theta_ref[1,1]", "theta_ref[1,2]") %in%
                    summ$variable))
  expect_true(all(c("sigma_y[1]", "sigma_y[2]") %in% summ$variable))
  expect_false(any(grepl("^theta_i\\[", summ$variable)))
  expect_false(any(grepl("^a_coef\\[", summ$variable)))

  out <- capture.output(print(fit))
  expect_true(any(grepl("p \\(theta_ref dim\\)\\s*:\\s*2", out)))
  expect_true(any(grepl("observations\\s*:\\s*300", out)))

  cf <- coef(fit)
  expect_s3_class(cf, "gdpar_coef")
  expect_equal(cf$p, 2L)
  expect_s3_class(cf$theta_ref, "data.frame")
  expect_named(cf$theta_ref, c("k", "mean", "q05", "q50", "q95"))
  expect_equal(cf$theta_ref$k, 1:2)
  expect_length(cf$a, 2L)
  for (k in 1:2) {
    expect_s3_class(cf$a[[k]], "data.frame")
    expect_named(cf$a[[k]],
                 c("term", "mean", "q05", "q50", "q95"))
    expect_equal(cf$a[[k]]$term, c("x1", "x2"))
  }
  expect_null(cf$b)
  expect_null(cf$W)
  cf_long <- as.data.frame(cf)
  expect_s3_class(cf_long, "data.frame")
  expect_named(cf_long, c("component", "g", "k", "identifier",
                          "x_name", "mean", "q05", "q50", "q95"))
  expect_equal(sum(cf_long$component == "theta_ref"), 2L)
  expect_equal(sum(cf_long$component == "a"), 4L)
})

test_that("gdpar() multi with W active: predict_from_newdata_multi reproduces Stan theta_i", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          "Set GDPAR_F6_FIT=1 to run the multi end-to-end smoke (compiles a Stan model; takes minutes).")
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  n <- 200L
  p <- 2L
  set.seed(7L)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  data_df <- data.frame(x1 = x1, x2 = x2)
  theta_ref_true <- c(0.4, -0.6)
  beta_a_true <- list(0.4, -0.3)
  sigma_true <- 0.3
  x1_c <- x1 - mean(x1)
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]] * x1_c,
    theta_ref_true[2] + beta_a_true[[2]] * x1_c
  )
  data_df$y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)
  spec <- amm_spec(
    p      = 2L,
    dims   = dimwise(a = ~ x1, b = NULL),
    W      = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1", "x2")
  )

  # spec has deliberate overlap (x1 in both a and x_vars). The
  # purpose of this test is the predict_from_newdata paridad, not
  # C4-bis. We skip the id-check explicitly; the post-fit
  # diagnostics still catch real non-identifiabilities via Rhat/ESS.
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family_multi("gaussian", p = 2L),
    amm           = spec,
    data          = data_df,
    chains        = 2L,
    iter_warmup   = 200L,
    iter_sampling = 200L,
    refresh       = 0L,
    seed          = 7L,
    skip_id_check = TRUE,
    verbose       = FALSE
  )

  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$p, 2L)
  expect_equal(fit$amm$W$dim, 4L)
  expect_equal(fit$amm$W$p, 2L)

  pred_draws <- predict(fit, summary = "draws")
  pred_newdata <- predict(fit, newdata = data_df, summary = "draws")
  expect_equal(dim(pred_draws), dim(pred_newdata))
  diff <- max(abs(pred_newdata - pred_draws))
  # Tolerance: 5e-6. The W block introduces FP reordering between
  # Stan's loop-based dot_product (C++) and R's crossprod (BLAS
  # vectorized). Empirically the diff lands in [1e-7, 2e-6] across
  # seeds and Stan template variants; 5e-6 gives a 2.5x safety
  # margin without losing fidelity (errors that matter would be
  # >1e-3).
  expect_lt(diff, 5e-6)

  # coef multi with W: confirms slot shape and per-k block structure.
  cf <- coef(fit)
  expect_s3_class(cf, "gdpar_coef")
  expect_length(cf$a, 2L)
  expect_length(cf$W, 2L)
  for (k in 1:2) {
    expect_s3_class(cf$a[[k]], "data.frame")
    expect_equal(cf$a[[k]]$term, "x1")
    expect_s3_class(cf$W[[k]], "data.frame")
    expect_named(cf$W[[k]],
                 c("basis_idx", "x_name", "mean", "q05", "q50", "q95"))
    expect_equal(sort(unique(cf$W[[k]]$basis_idx)), 1:2)
    expect_equal(sort(unique(cf$W[[k]]$x_name)), c("x1", "x2"))
  }
})

test_that("gdpar() multi auto-promotion end-to-end: univariate family becomes multi inside fit", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          "Set GDPAR_F6_FIT=1 to run the multi end-to-end smoke (compiles a Stan model; takes minutes).")
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  n <- 150L
  p <- 2L
  set.seed(11L)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  data_df <- data.frame(x1 = x1, x2 = x2)
  theta_ref_true <- c(0.5, -0.5)
  beta_a_true <- list(c(0.5, -0.5), c(-0.4, 0.4))
  sigma_true <- 0.3
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1_c + beta_a_true[[1]][2] * x2_c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1_c + beta_a_true[[2]][2] * x2_c
  )
  data_df$y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))

  suppressWarnings(
    expect_message(
      fit <- gdpar(
        formula       = y ~ x1 + x2,
        family        = gdpar_family("gaussian"),
        amm           = spec,
        data          = data_df,
        chains        = 2L,
        iter_warmup   = 200L,
        iter_sampling = 200L,
        refresh       = 0L,
        seed          = 11L,
        verbose       = TRUE
      ),
      regexp = "Auto-promoted univariate family 'gaussian'"
    )
  )
  expect_s3_class(fit$family, "gdpar_family_multi")
  expect_equal(fit$family$p, 2L)
})

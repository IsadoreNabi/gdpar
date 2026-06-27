# Sub-phase 8.3.9 (2026-05-22): re-validation of the S3 methods
# predict.gdpar_fit / summary.gdpar_fit / coef.gdpar_fit /
# gdpar_loo() against fits with K > 1, p > 1 and bspline W bases
# (D5 = 5D: smokes + tests + fuzz, closes gotcha 48).
#
# Structural tests (no Stan): use mock fit objects with synthetic
# draws to exercise dispatch and shape contracts. Gated smokes
# (GDPAR_RUN_STAN_S3_FUZZ) exercise the real path on a small fit.

# ---- Mock helpers -------------------------------------------------------

.make_K_mock_fit <- function(K = 2L, n = 20L, S = 100L, seed = 99L,
                             family_name = "gaussian") {
  set.seed(seed)
  # Build draws_matrix-like object containing theta_ref, theta_i_k,
  # eta, log_lik, y_pred per the K-individual Stan template.
  var_names <- c()
  col_data  <- list()
  # theta_ref[k]
  for (k in seq_len(K)) {
    var_names <- c(var_names, sprintf("theta_ref[%d]", k))
    col_data[[length(col_data) + 1L]] <- rnorm(S)
  }
  # theta_i_k[i, k]
  for (k in seq_len(K)) {
    for (i in seq_len(n)) {
      var_names <- c(var_names, sprintf("theta_i_k[%d,%d]", i, k))
      col_data[[length(col_data) + 1L]] <- rnorm(S)
    }
  }
  # eta_k[i, k]
  for (k in seq_len(K)) {
    for (i in seq_len(n)) {
      var_names <- c(var_names, sprintf("eta_k[%d,%d]", i, k))
      col_data[[length(col_data) + 1L]] <- rnorm(S)
    }
  }
  # log_lik[i]
  for (i in seq_len(n)) {
    var_names <- c(var_names, sprintf("log_lik[%d]", i))
    col_data[[length(col_data) + 1L]] <- rnorm(S)
  }
  # y_pred[i]
  for (i in seq_len(n)) {
    var_names <- c(var_names, sprintf("y_pred[%d]", i))
    col_data[[length(col_data) + 1L]] <- rnorm(S)
  }
  mat <- do.call(cbind, col_data)
  colnames(mat) <- var_names
  draws_arr <- posterior::as_draws_matrix(mat)
  fake_fit <- list()
  fake_fit$draws <- function(variables = NULL, format = "draws_matrix") {
    if (is.null(variables)) {
      return(draws_arr)
    }
    keep <- character(0)
    for (v in variables) {
      keep <- c(keep, var_names[startsWith(var_names, paste0(v, "["))],
                if (v %in% var_names) v else character(0))
    }
    keep <- unique(keep)
    posterior::subset_draws(draws_arr, variable = keep)
  }
  fam <- gdpar_family(family_name)
  obj <- list(
    fit                = fake_fit,
    amm_list_canonical = list(),
    family             = fam,
    prior              = gdpar_prior(),
    design_K           = list(),
    anchor             = rep(0, K),
    stan_data          = list(y_real = rnorm(n), y_int = integer(n)),
    K                  = K,
    slot_names         = if (K == 2L) c("mu", "sigma")
                         else paste0("slot_", seq_len(K)),
    diagnostics        = list(converged = TRUE,
                              rhat_max = 1.001,
                              ess_bulk_min = 500,
                              divergent_count = 0L),
    path               = "bayes"
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

# ---- Section 1: predict.gdpar_fit on K > 1 (Ruta B 8.3.9) ----------------

test_that("predict.gdpar_fit on K>1 mock returns array S*n*K (default type='theta_i', summary='draws')", {
  obj <- .make_K_mock_fit(K = 2L, n = 10L, S = 50L)
  arr <- predict(obj)
  expect_true(is.array(arr))
  expect_equal(dim(arr), c(50L, 10L, 2L))
  expect_identical(dimnames(arr)[[3L]], c("mu", "sigma"))
})

test_that("predict.gdpar_fit on K=3 mock returns array S*n*3", {
  obj <- .make_K_mock_fit(K = 3L, n = 15L, S = 40L,
                          family_name = "student_t")
  arr <- predict(obj)
  expect_equal(dim(arr), c(40L, 15L, 3L))
})

test_that("predict.gdpar_fit on K>1 supports type='response' applying per-slot inv_link", {
  obj <- .make_K_mock_fit(K = 2L, n = 8L, S = 30L)
  arr_lp <- predict(obj, type = "linear_predictor")
  arr_resp <- predict(obj, type = "response")
  # For Gaussian K=2 (location identity, sigma log), inv_link of slot 1
  # is identity, inv_link of slot 2 is exp; so arr_resp[, , 1] equals
  # arr_lp[, , 1] but arr_resp[, , 2] = exp(arr_lp[, , 2]).
  expect_equal(arr_resp[, , 1L], arr_lp[, , 1L])
  expect_equal(arr_resp[, , 2L], exp(arr_lp[, , 2L]))
})

test_that("predict.gdpar_fit on K>1 with summary='mean_se' returns named list of length K", {
  obj <- .make_K_mock_fit(K = 2L, n = 12L, S = 40L)
  out <- predict(obj, summary = "mean_se")
  expect_true(is.list(out))
  expect_length(out, 2L)
  expect_identical(names(out), c("mu", "sigma"))
  for (k in seq_len(2L)) {
    expect_true(is.data.frame(out[[k]]))
    expect_equal(nrow(out[[k]]), 12L)
    expect_identical(names(out[[k]]), c("mean", "se"))
  }
})

test_that("predict.gdpar_fit on K>1 with summary='quantiles' returns named list with q05/q50/q95", {
  obj <- .make_K_mock_fit(K = 2L, n = 10L, S = 40L)
  out <- predict(obj, summary = "quantiles")
  expect_length(out, 2L)
  for (k in seq_len(2L)) {
    expect_identical(names(out[[k]]), c("q05", "q50", "q95"))
  }
})

# predict.gdpar_fit K>1 with newdata is implemented in 8.3.10 via
# predict_from_newdata_K. Two unsupported branches remain (queued for
# Session 8.4): K-individual with grouping (J_groups>1), and K-individual
# with B-spline W bases on newdata. Polynomial W and ungrouped K>1 are
# supported.

# Mock helper used by the predict-newdata K>1 tests: populates
# theta_ref_k[1,k], a_coef_k[k,j], c_b_k[k,j], W_raw[r,c], sigma_W[1]
# plus a *minimal* design_K (means + sds + per-slot column-name lists)
# and amm_list_canonical with explicit ~ x1 / ~ x2 formulas so that the
# per-slot design rebuild via stats::model.matrix is exercised. Anchor
# is set to 0 for every slot.
.make_K_mock_fit_predict_newdata <- function(K = 2L, S = 80L,
                                              seed = 717L,
                                              a_per_k = c(1L, 1L),
                                              W_dim = 0L,
                                              W_x = character(0),
                                              cp_W = TRUE) {
  set.seed(seed)
  stopifnot(length(a_per_k) == K)
  J_a_max <- max(a_per_k, 1L)
  d_x <- length(W_x)
  var_names <- c()
  col_data <- list()
  for (k in seq_len(K)) {
    var_names <- c(var_names, sprintf("theta_ref_k[%d,%d]", 1L, k))
    col_data[[length(col_data) + 1L]] <- rnorm(S, mean = 0.1 * k, sd = 0.05)
  }
  for (k in seq_len(K)) {
    for (j in seq_len(J_a_max)) {
      var_names <- c(var_names, sprintf("a_coef_k[%d,%d]", k, j))
      if (j <= a_per_k[k]) {
        col_data[[length(col_data) + 1L]] <- rnorm(S, mean = 0.5, sd = 0.05)
      } else {
        col_data[[length(col_data) + 1L]] <- rep(0, S)
      }
    }
  }
  for (k in seq_len(K)) {
    var_names <- c(var_names, sprintf("c_b_k[%d,%d]", k, 1L))
    col_data[[length(col_data) + 1L]] <- rep(0, S)
  }
  if (W_dim > 0L && d_x > 0L) {
    for (r in seq_len(W_dim)) {
      for (c_idx in seq_len(d_x)) {
        var_names <- c(var_names, sprintf("W_raw[%d,%d]", r, c_idx))
        col_data[[length(col_data) + 1L]] <- rnorm(S, sd = 0.2)
      }
    }
    if (!cp_W) {
      var_names <- c(var_names, "sigma_W[1]")
      col_data[[length(col_data) + 1L]] <- abs(rnorm(S, 0.5, 0.1))
    }
  }
  mat <- do.call(cbind, col_data)
  colnames(mat) <- var_names
  draws_arr <- posterior::as_draws_matrix(mat)
  fake_fit <- list()
  fake_fit$draws <- function(variables = NULL, format = "draws_matrix") {
    if (is.null(variables)) return(draws_arr)
    keep <- character(0)
    for (v in variables) {
      keep <- c(keep, var_names[startsWith(var_names, paste0(v, "["))],
                if (v %in% var_names) v else character(0))
    }
    keep <- unique(keep)
    posterior::subset_draws(draws_arr, variable = keep)
  }
  slot_names <- if (K == 2L) c("mu", "sigma") else
    paste0("slot_", seq_len(K))
  Z_a_k_means_list <- lapply(seq_len(K), function(k) {
    if (a_per_k[k] > 0L) rep(0, a_per_k[k]) else numeric(0)
  })
  Z_a_k_names_list <- lapply(seq_len(K), function(k) {
    if (a_per_k[k] > 0L) sprintf("z_a_%s_%d", slot_names[k],
                                 seq_len(a_per_k[k]))
    else character(0)
  })
  amm_list_canonical <- lapply(seq_len(K), function(k) {
    a_form <- if (a_per_k[k] > 0L)
      stats::as.formula(paste0("~ x", k)) else NULL
    W_spec <- if (W_dim > 0L && d_x > 0L && k == 1L)
      list(type = "polynomial", dim = as.integer(W_dim)) else NULL
    out <- list(a = a_form, b = NULL, W = W_spec)
    class(out) <- c("amm_spec", "list")
    out
  })
  names(amm_list_canonical) <- slot_names
  design_K <- list(
    K = K, slot_names = slot_names,
    Z_a_k_names_list = Z_a_k_names_list,
    Z_b_k_names_list = lapply(seq_len(K), function(k) character(0)),
    Z_a_k_means_list = Z_a_k_means_list,
    Z_b_k_means_list = lapply(seq_len(K), function(k) numeric(0)),
    X_names = W_x,
    X_means = if (d_x > 0L) rep(0, d_x) else numeric(0),
    X_sds   = if (d_x > 0L) rep(1, d_x) else numeric(0)
  )
  fam <- gdpar_family("gaussian")
  obj <- list(
    fit                = fake_fit,
    amm_list_canonical = amm_list_canonical,
    family             = fam,
    prior              = gdpar_prior(),
    design_K           = design_K,
    anchor             = rep(0, K),
    stan_data          = list(),
    parametrization    = list(cp_W = cp_W),
    K                  = K,
    slot_names         = slot_names,
    diagnostics        = list(converged = TRUE, rhat_max = 1.001,
                              ess_bulk_min = 500, divergent_count = 0L),
    path               = "bayes"
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

test_that("predict.gdpar_fit on K>1 with newdata returns array S*n_new*K (polynomial W=0)", {
  obj <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                           a_per_k = c(1L, 1L))
  nd <- data.frame(x1 = rnorm(8L), x2 = rnorm(8L))
  arr <- predict(obj, newdata = nd)
  expect_true(is.array(arr))
  expect_equal(dim(arr), c(60L, 8L, 2L))
  expect_identical(dimnames(arr)[[3L]], c("mu", "sigma"))
})

test_that("predict.gdpar_fit on K>1 with newdata + W (polynomial) injects modulating contribution", {
  obj <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                           a_per_k = c(1L, 0L),
                                           W_dim = 2L,
                                           W_x = c("xw"),
                                           cp_W = TRUE)
  nd <- data.frame(x1 = rnorm(6L), xw = rnorm(6L))
  arr <- predict(obj, newdata = nd)
  expect_equal(dim(arr), c(60L, 6L, 2L))
  # W is globally shared so slot sigma (no slot-declared W) also
  # receives the modulating contribution; eta_sigma is not just
  # theta_ref_k[, sigma].
  obj2 <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                            a_per_k = c(1L, 0L),
                                            W_dim = 0L,
                                            W_x = character(0),
                                            cp_W = TRUE, seed = 717L)
  arr2 <- predict(obj2, newdata = nd[, "x1", drop = FALSE])
  # slot sigma differs between obj (with W) and obj2 (without W)
  expect_false(isTRUE(all.equal(arr[, , 2L], arr2[, , 2L])))
})

test_that("predict.gdpar_fit on K>1 with newdata + cp_W=FALSE multiplies by sigma_W", {
  obj1 <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                            a_per_k = c(1L, 0L),
                                            W_dim = 1L, W_x = c("xw"),
                                            cp_W = TRUE, seed = 444L)
  obj2 <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                            a_per_k = c(1L, 0L),
                                            W_dim = 1L, W_x = c("xw"),
                                            cp_W = FALSE, seed = 444L)
  nd <- data.frame(x1 = rnorm(5L), xw = rnorm(5L))
  arr1 <- predict(obj1, newdata = nd)
  arr2 <- predict(obj2, newdata = nd)
  expect_false(isTRUE(all.equal(arr1, arr2)))
})

test_that("predict.gdpar_fit on K>1 with newdata supports summary='mean_se' per slot", {
  obj <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                           a_per_k = c(1L, 1L))
  nd <- data.frame(x1 = rnorm(7L), x2 = rnorm(7L))
  out <- predict(obj, newdata = nd, summary = "mean_se")
  expect_true(is.list(out))
  expect_length(out, 2L)
  expect_identical(names(out), c("mu", "sigma"))
  for (k in seq_along(out)) {
    expect_true(is.data.frame(out[[k]]))
    expect_equal(nrow(out[[k]]), 7L)
    expect_identical(names(out[[k]]), c("mean", "se"))
  }
})

test_that("predict.gdpar_fit on K>1 with newdata + grouping raises gdpar_unsupported_feature_error", {
  obj <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                           a_per_k = c(1L, 1L))
  obj$group_info <- list(levels = c("A", "B"))
  nd <- data.frame(x1 = rnorm(4L), x2 = rnorm(4L))
  err <- tryCatch(predict(obj, newdata = nd), error = function(e) e)
  expect_true(inherits(err, "gdpar_unsupported_feature_error"))
  expect_match(conditionMessage(err), "grouping|J_groups",
               ignore.case = TRUE)
})

test_that("predict.gdpar_fit on K>1 with newdata + bspline W raises gdpar_unsupported_feature_error", {
  obj <- .make_K_mock_fit_predict_newdata(K = 2L, S = 60L,
                                           a_per_k = c(1L, 0L),
                                           W_dim = 3L, W_x = c("xw"))
  # Swap W spec to bspline post-hoc to exercise the guard
  obj$amm_list_canonical[[1L]]$W$type <- "bspline"
  nd <- data.frame(x1 = rnorm(4L), xw = rnorm(4L))
  err <- tryCatch(predict(obj, newdata = nd), error = function(e) e)
  expect_true(inherits(err, "gdpar_unsupported_feature_error"))
  expect_match(conditionMessage(err), "B-spline|bspline",
               ignore.case = TRUE)
})

# ---- Section 2: summary.gdpar_fit on K > 1 ------------------------------

test_that("summary.gdpar_fit on K>1 mock returns a data.frame with theta_ref entries", {
  skip_if_not_installed("posterior")
  obj <- .make_K_mock_fit(K = 2L, n = 15L, S = 80L)
  s <- summary(obj)
  expect_true(is.data.frame(s) || inherits(s, "tbl_df"))
  var_col <- if ("variable" %in% names(s)) s$variable else rownames(s)
  expect_true(any(grepl("^theta_ref", var_col)))
})

test_that("summary.gdpar_fit excludes auxiliary draws (eta/log_lik/y_pred/theta_i)", {
  skip_if_not_installed("posterior")
  obj <- .make_K_mock_fit(K = 2L, n = 12L, S = 60L)
  s <- summary(obj)
  var_col <- if ("variable" %in% names(s)) s$variable else rownames(s)
  expect_false(any(grepl("^eta_k", var_col)))
  expect_false(any(grepl("^log_lik", var_col)))
  expect_false(any(grepl("^y_pred", var_col)))
  expect_false(any(grepl("^theta_i_k", var_col)))
})

# ---- Section 2.5: coef.gdpar_fit on K > 1 (decision E4.A 8.3.10) --------

# Extended mock helper for coef.gdpar_fit K>1 tests. Populates the
# K-individual Stan variables that coef() needs: theta_ref_k[1, k] (one
# group, K slots), a_coef_k[k, j] padded to max_J_a, c_b_k[k, j] padded
# to max_J_b, W_raw[r, c] and sigma_W[1]. Also wires design_K and
# amm_list_canonical so that per-slot column-name and AMM-component
# extraction works.
.make_K_mock_fit_for_coef <- function(K = 2L, n = 20L, S = 80L,
                                       seed = 313L,
                                       a_per_k = c(2L, 1L),
                                       b_per_k = c(0L, 0L),
                                       W_dim = 0L, W_x = character(0),
                                       cp_W = TRUE) {
  set.seed(seed)
  stopifnot(length(a_per_k) == K, length(b_per_k) == K)
  J_a_max <- max(a_per_k, 1L)
  J_b_max <- max(b_per_k, 1L)
  d_x <- length(W_x)
  var_names <- c()
  col_data  <- list()
  # theta_ref_k[1, k]: array[J_groups=1] vector[K]
  for (k in seq_len(K)) {
    var_names <- c(var_names, sprintf("theta_ref_k[%d,%d]", 1L, k))
    col_data[[length(col_data) + 1L]] <- rnorm(S, mean = 0.1 * k, sd = 0.05)
  }
  # a_coef_k[k, j]: padded to J_a_max; padding entries are zeros
  for (k in seq_len(K)) {
    for (j in seq_len(J_a_max)) {
      var_names <- c(var_names, sprintf("a_coef_k[%d,%d]", k, j))
      if (j <= a_per_k[k]) {
        col_data[[length(col_data) + 1L]] <- rnorm(S, mean = 0.4, sd = 0.1)
      } else {
        col_data[[length(col_data) + 1L]] <- rep(0, S)
      }
    }
  }
  # c_b_k[k, j]: padded to J_b_max
  for (k in seq_len(K)) {
    for (j in seq_len(J_b_max)) {
      var_names <- c(var_names, sprintf("c_b_k[%d,%d]", k, j))
      if (j <= b_per_k[k]) {
        col_data[[length(col_data) + 1L]] <- rnorm(S, mean = 0.2, sd = 0.05)
      } else {
        col_data[[length(col_data) + 1L]] <- rep(0, S)
      }
    }
  }
  # W_raw[r, c] + sigma_W[1] (global; only if W_dim > 0)
  if (W_dim > 0L && d_x > 0L) {
    for (r in seq_len(W_dim)) {
      for (c_idx in seq_len(d_x)) {
        var_names <- c(var_names, sprintf("W_raw[%d,%d]", r, c_idx))
        col_data[[length(col_data) + 1L]] <- rnorm(S, sd = 0.3)
      }
    }
    if (!cp_W) {
      var_names <- c(var_names, "sigma_W[1]")
      col_data[[length(col_data) + 1L]] <- abs(rnorm(S, mean = 0.5, sd = 0.1))
    }
  }
  mat <- do.call(cbind, col_data)
  colnames(mat) <- var_names
  draws_arr <- posterior::as_draws_matrix(mat)
  fake_fit <- list()
  fake_fit$draws <- function(variables = NULL, format = "draws_matrix") {
    if (is.null(variables)) return(draws_arr)
    keep <- character(0)
    for (v in variables) {
      keep <- c(keep, var_names[startsWith(var_names, paste0(v, "["))],
                if (v %in% var_names) v else character(0))
    }
    keep <- unique(keep)
    posterior::subset_draws(draws_arr, variable = keep)
  }
  slot_names <- if (K == 2L) c("mu", "sigma") else
    paste0("slot_", seq_len(K))
  Z_a_k_names_list <- lapply(seq_len(K), function(k) {
    if (a_per_k[k] > 0L) sprintf("z_a_%s_%d", slot_names[k],
                                 seq_len(a_per_k[k]))
    else character(0)
  })
  Z_b_k_names_list <- lapply(seq_len(K), function(k) {
    if (b_per_k[k] > 0L) sprintf("z_b_%s_%d", slot_names[k],
                                 seq_len(b_per_k[k]))
    else character(0)
  })
  amm_list_canonical <- lapply(seq_len(K), function(k) {
    out <- list(
      a = if (a_per_k[k] > 0L) stats::as.formula("~ x1") else NULL,
      b = if (b_per_k[k] > 0L) stats::as.formula("~ x2") else NULL,
      W = if (W_dim > 0L && d_x > 0L && k == 1L)
        list(type = "polynomial", dim = as.integer(W_dim)) else NULL
    )
    class(out) <- c("amm_spec", "list")
    out
  })
  names(amm_list_canonical) <- slot_names
  design_K <- list(
    K = K, slot_names = slot_names,
    Z_a_k_names_list = Z_a_k_names_list,
    Z_b_k_names_list = Z_b_k_names_list,
    X_names = W_x
  )
  fam <- if (K == 2L) gdpar_family("gaussian") else gdpar_family("student_t")
  obj <- list(
    fit                = fake_fit,
    amm_list_canonical = amm_list_canonical,
    family             = fam,
    prior              = gdpar_prior(),
    design_K           = design_K,
    anchor             = rep(0, K),
    stan_data          = list(y_real = rnorm(n)),
    parametrization    = list(cp_W = cp_W),
    K                  = K,
    slot_names         = slot_names,
    diagnostics        = list(converged = TRUE, rhat_max = 1.001,
                              ess_bulk_min = 500, divergent_count = 0L),
    path               = "bayes"
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

test_that("coef.gdpar_fit on K>1 mock returns named list of length K", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                    a_per_k = c(2L, 1L))
  co <- coef(obj)
  expect_true(is.list(co))
  expect_length(co, 2L)
  expect_identical(names(co), c("mu", "sigma"))
  for (k in seq_along(co)) {
    expect_s3_class(co[[k]], "gdpar_coef")
    expect_equal(co[[k]]$p, 1L)
    expect_equal(nrow(co[[k]]$theta_ref), 1L)
  }
})

test_that("coef.gdpar_fit on K>1 reports per-slot a coefficients with the expected term names", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                    a_per_k = c(2L, 1L))
  co <- coef(obj)
  expect_true(is.list(co$mu$a))
  expect_length(co$mu$a, 1L)
  expect_identical(co$mu$a[[1L]]$term, c("z_a_mu_1", "z_a_mu_2"))
  expect_true(is.list(co$sigma$a))
  expect_length(co$sigma$a, 1L)
  expect_identical(co$sigma$a[[1L]]$term, "z_a_sigma_1")
})

test_that("coef.gdpar_fit on K>1 marks NULL for slots without an active component", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                    a_per_k = c(1L, 0L),
                                    b_per_k = c(0L, 0L))
  co <- coef(obj)
  expect_null(co$sigma$a)
  expect_null(co$mu$b)
  expect_null(co$sigma$b)
  expect_null(co$mu$W)
  expect_null(co$sigma$W)
})

test_that("coef.gdpar_fit on K>1 reports b coefficients via c_b_k draws", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                    a_per_k = c(0L, 0L),
                                    b_per_k = c(1L, 1L))
  co <- coef(obj)
  expect_true(is.list(co$mu$b))
  expect_identical(co$mu$b[[1L]]$term, "z_b_mu_1")
  expect_identical(co$sigma$b[[1L]]$term, "z_b_sigma_1")
})

test_that("coef.gdpar_fit on K>1 with W_global attaches to every slot (canonical decision: W is globally shared)", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                    a_per_k = c(1L, 0L),
                                    W_dim = 2L, W_x = c("x1", "x2"),
                                    cp_W = TRUE)
  co <- coef(obj)
  for (slot in c("mu", "sigma")) {
    expect_true(is.list(co[[slot]]$W),
                info = paste("slot =", slot))
    expect_length(co[[slot]]$W, 1L)
    expect_true(is.data.frame(co[[slot]]$W[[1L]]))
    expect_identical(sort(unique(co[[slot]]$W[[1L]]$x_name)),
                     c("x1", "x2"))
    expect_equal(nrow(co[[slot]]$W[[1L]]), 2L * 2L)
  }
  # The W block is the same across slots (bit-exact copy of the
  # global posterior); means must coincide.
  expect_equal(co$mu$W[[1L]]$mean, co$sigma$W[[1L]]$mean)
})

test_that("coef.gdpar_fit on K>1 with W and cp_W=FALSE multiplies by sigma_W", {
  obj1 <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                     a_per_k = c(1L, 0L),
                                     W_dim = 2L, W_x = c("x1"),
                                     cp_W = TRUE, seed = 999L)
  obj2 <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                     a_per_k = c(1L, 0L),
                                     W_dim = 2L, W_x = c("x1"),
                                     cp_W = FALSE, seed = 999L)
  co1 <- coef(obj1)
  co2 <- coef(obj2)
  # cp_W=FALSE multiplies by sigma_W per draw; means differ across paths.
  expect_false(isTRUE(all.equal(co1$mu$W[[1L]]$mean, co2$mu$W[[1L]]$mean)))
})

test_that("coef.gdpar_fit on K>1 with grouping raises gdpar_unsupported_feature_error (deuda Session 8.4)", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L)
  obj$group_info <- list(levels = c("A", "B"))
  err <- tryCatch(coef(obj), error = function(e) e)
  expect_true(inherits(err, "gdpar_unsupported_feature_error"))
  expect_match(conditionMessage(err), "grouping|J_groups",
               ignore.case = TRUE)
})

test_that("coef.gdpar_fit on K>1 supports as.data.frame per slot", {
  obj <- .make_K_mock_fit_for_coef(K = 2L, n = 20L, S = 80L,
                                    a_per_k = c(1L, 1L))
  co <- coef(obj)
  df_list <- lapply(co, as.data.frame)
  expect_true(all(vapply(df_list, is.data.frame, logical(1L))))
  expect_true(all(vapply(df_list, function(d) nrow(d) >= 2L, logical(1L))))
})

# ---- Section 3: print.gdpar_fit on K > 1 --------------------------------

test_that("print.gdpar_fit on K>1 mock produces visible output without erroring", {
  obj <- .make_K_mock_fit(K = 3L, n = 10L, S = 30L)
  expect_output(print(obj), "gdpar_fit", fixed = TRUE)
})

# ---- Section 4: gdpar_loo on K > 1 --------------------------------------
# Note: gdpar_loo reads object$amm$p (NULL for K-individual fits, so
# treated as p = 1). log_lik in amm_distrib_K.stan is vector[n], so the
# scalar aggregation path is the correct one.

test_that("gdpar_loo on K>1 mock matches the p=1 aggregation contract", {
  skip_if_not_installed("loo")
  obj <- .make_K_mock_fit(K = 2L, n = 10L, S = 40L)
  # gdpar_loo expects fit$fit$draws(variables = "log_lik", format = "draws_array")
  obj$fit$draws <- function(variables = NULL, format = "draws_array") {
    var_names <- paste0("log_lik[", seq_len(10L), "]")
    mat <- matrix(rnorm(40L * 10L), nrow = 40L)
    colnames(mat) <- var_names
    draws_arr <- posterior::as_draws_array(
      posterior::as_draws_matrix(mat)
    )
    if (!is.null(variables) && "log_lik" %in% variables) {
      return(draws_arr)
    }
    draws_arr
  }
  lo <- tryCatch(suppressWarnings(gdpar_loo(obj)), error = function(e) e)
  expect_false(inherits(lo, "error"))
  expect_true(inherits(lo, "psis_loo") || inherits(lo, "loo"))
})

# ---- Section 5: fuzz multiconf (gated) ----------------------------------
# Sets GDPAR_RUN_STAN_S3_FUZZ=1 to compile and fit several small
# configurations and exercise predict / summary / coef / gdpar_loo
# against actual posterior draws. D5 = 5D ranqueada por maxima
# robustez en 8.3.9 scoping.

.skip_if_no_s3_fuzz <- function() {
  if (!nzchar(Sys.getenv("GDPAR_RUN_STAN_S3_FUZZ"))) {
    skip("GDPAR_RUN_STAN_S3_FUZZ not set; skipping S3 methods fuzz")
  }
  skip_if_not_installed("cmdstanr")
}

.fuzz_K1_gaussian <- function() {
  set.seed(101L)
  n <- 60L
  x <- rnorm(n); z <- rnorm(n)
  y <- 0.4 + 0.7 * (x - mean(x)) + rnorm(n, 0, 0.5)
  d <- data.frame(y = y, x = x, z = z)
  suppressMessages(suppressWarnings(
    gdpar(y ~ a(x), data = d, family = gdpar_family("gaussian"),
          chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
          refresh = 0L)
  ))
}

.fuzz_K2_gaussian <- function() {
  set.seed(202L)
  n <- 80L
  x1 <- rnorm(n); x2 <- rnorm(n)
  mu_true   <- 0.5 + 0.7 * (x1 - mean(x1))
  log_sigma <- -0.2 + 0.5 * (x2 - mean(x2))
  y <- rnorm(n, mu_true, exp(log_sigma))
  d <- data.frame(y = y, x1 = x1, x2 = x2)
  fs <- gdpar_bf(y ~ a(x1), sigma ~ a(x2))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("gaussian"),
          chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
          refresh = 0L)
  ))
}

test_that("FUZZ: predict / summary / coef / gdpar_loo on K=1+p=1 Gaussian", {
  .skip_if_no_s3_fuzz()
  fit <- .fuzz_K1_gaussian()
  expect_silent(s <- summary(fit))
  expect_true(is.data.frame(s) || inherits(s, "tbl_df"))
  pred_draws <- predict(fit, type = "theta_i", summary = "draws")
  expect_true(is.matrix(pred_draws))
  expect_equal(ncol(pred_draws), 60L)
  co <- tryCatch(coef(fit), error = function(e) e)
  expect_false(inherits(co, "error"))
  lo <- tryCatch(gdpar_loo(fit), error = function(e) e)
  expect_false(inherits(lo, "error"))
  expect_true(inherits(lo, "psis_loo") || inherits(lo, "loo"))
})

test_that("FUZZ: predict / summary / coef / gdpar_loo on K=2 Gaussian dist regression", {
  .skip_if_no_s3_fuzz()
  fit <- .fuzz_K2_gaussian()
  expect_silent(s <- summary(fit))
  expect_true(is.data.frame(s) || inherits(s, "tbl_df"))
  # predict on K>1 is part of gotcha 48 (still open in 8.3.9): it
  # should either dispatch correctly or error informatively.
  pred_or_err <- tryCatch(predict(fit), error = function(e) e)
  if (inherits(pred_or_err, "error")) {
    expect_match(conditionMessage(pred_or_err),
                 "theta_i|K|coordinate|individual",
                 ignore.case = TRUE)
  } else {
    expect_true(is.matrix(pred_or_err) || is.array(pred_or_err))
  }
  co <- tryCatch(coef(fit), error = function(e) e)
  expect_false(inherits(co, "error"))
  lo <- tryCatch(gdpar_loo(fit), error = function(e) e)
  expect_false(inherits(lo, "error"))
  expect_true(inherits(lo, "psis_loo") || inherits(lo, "loo"))
})

# ---- Section 6: residuals.gdpar_fit on K>1 fuzz -------------------------

test_that("FUZZ: residuals.gdpar_fit on K=2 Gaussian dist regression", {
  .skip_if_no_s3_fuzz()
  fit <- .fuzz_K2_gaussian()
  for (t in c("response", "pearson", "deviance", "quantile")) {
    r <- residuals(fit, type = t, randomize_seed = 1L)
    expect_length(r, 80L)
    expect_true(all(is.finite(r)),
                info = sprintf("type=%s", t))
  }
})

test_that("FUZZ: gdpar_posterior_predict on K=2 Gaussian dist regression", {
  .skip_if_no_s3_fuzz()
  fit <- .fuzz_K2_gaussian()
  pp <- gdpar_posterior_predict(fit)
  expect_true(is.matrix(pp))
  expect_equal(ncol(pp), 80L)
})

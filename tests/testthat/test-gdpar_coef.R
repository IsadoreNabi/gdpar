make_theta_df <- function(p = 2L, seed = 1L) {
  set.seed(seed)
  S <- 60L
  dr <- matrix(rnorm(S * p, mean = seq_len(p)), nrow = S, ncol = p,
               byrow = TRUE)
  build_coef_theta_ref_df(dr, p)
}

make_term_df <- function(n_terms = 2L, seed = 11L) {
  set.seed(seed)
  S <- 60L
  dr <- matrix(rnorm(S * n_terms), nrow = S, ncol = n_terms)
  build_coef_term_df(dr, paste0("term", seq_len(n_terms)))
}

make_W_df <- function(basis_dim = 2L, d_x = 2L, seed = 21L) {
  set.seed(seed)
  S <- 60L
  dr <- matrix(rnorm(S * basis_dim * d_x), nrow = S,
               ncol = basis_dim * d_x)
  build_coef_W_df(dr, basis_dim, paste0("x", seq_len(d_x)))
}

# ---- helpers --------------------------------------------------------

test_that("build_coef_theta_ref_df produces the expected structure", {
  set.seed(1)
  dr <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  out <- build_coef_theta_ref_df(dr, p = 3L)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("k", "mean", "q05", "q50", "q95"))
  expect_equal(out$k, 1:3)
  expect_equal(nrow(out), 3L)
  expect_true(all(is.finite(out$mean)))
})

test_that("build_coef_theta_ref_df aborts when draws cols != p", {
  dr <- matrix(rnorm(50 * 2), nrow = 50, ncol = 2)
  expect_error(
    build_coef_theta_ref_df(dr, p = 3L),
    class = "gdpar_internal_error"
  )
})

test_that("build_coef_term_df returns tidy term-summary df", {
  out <- make_term_df(n_terms = 3L)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "mean", "q05", "q50", "q95"))
  expect_equal(out$term, c("term1", "term2", "term3"))
  expect_equal(nrow(out), 3L)
})

test_that("build_coef_term_df returns NULL on empty draws", {
  expect_null(build_coef_term_df(NULL, character(0)))
  expect_null(build_coef_term_df(matrix(numeric(0), nrow = 0, ncol = 0),
                                 character(0)))
})

test_that("build_coef_term_df aborts when term_names mismatch cols", {
  dr <- matrix(rnorm(50 * 2), nrow = 50, ncol = 2)
  expect_error(
    build_coef_term_df(dr, term_names = c("a", "b", "c")),
    class = "gdpar_internal_error"
  )
})

test_that("build_coef_W_df returns tidy W-summary df with row-major flat layout", {
  out <- make_W_df(basis_dim = 3L, d_x = 2L)
  expect_s3_class(out, "data.frame")
  expect_named(out,
               c("basis_idx", "x_name", "mean", "q05", "q50", "q95"))
  expect_equal(nrow(out), 6L)
  expect_equal(out$basis_idx, rep(1:3, times = 2))
  expect_equal(out$x_name, rep(c("x1", "x2"), each = 3))
})

test_that("build_coef_W_df returns NULL on empty inputs", {
  expect_null(build_coef_W_df(NULL, 2L, c("x1")))
  expect_null(build_coef_W_df(matrix(0, nrow = 5L, ncol = 0L),
                              basis_dim = 0L, x_names = character(0)))
})

test_that("build_coef_W_df aborts on shape mismatch", {
  dr <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  expect_error(
    build_coef_W_df(dr, basis_dim = 3L, x_names = c("x1", "x2")),
    class = "gdpar_internal_error"
  )
})

# ---- constructor + validators --------------------------------------

test_that("new_gdpar_coef builds a minimal valid object (p = 1, all components)", {
  theta_df <- make_theta_df(p = 1L)
  a_df <- make_term_df(n_terms = 2L, seed = 31L)
  b_df <- make_term_df(n_terms = 1L, seed = 32L)
  W_df <- make_W_df(basis_dim = 2L, d_x = 2L, seed = 33L)
  obj <- new_gdpar_coef(
    theta_ref = theta_df,
    a = list(a_df), b = list(b_df), W = list(W_df), p = 1L
  )
  expect_s3_class(obj, "gdpar_coef")
  expect_equal(obj$p, 1L)
  expect_equal(obj$summary_stats, c("mean", "q05", "q50", "q95"))
})

test_that("new_gdpar_coef accepts NULL components (absent at AMM level)", {
  theta_df <- make_theta_df(p = 2L)
  obj <- new_gdpar_coef(theta_ref = theta_df, p = 2L)
  expect_s3_class(obj, "gdpar_coef")
  expect_null(obj$a); expect_null(obj$b); expect_null(obj$W)
})

test_that("new_gdpar_coef accepts NULL per-coord entries (coord inactive)", {
  theta_df <- make_theta_df(p = 2L)
  a_df <- make_term_df(seed = 41L)
  obj <- new_gdpar_coef(
    theta_ref = theta_df,
    a = list(a_df, NULL), p = 2L
  )
  expect_null(obj$a[[2]])
  expect_s3_class(obj$a[[1]], "data.frame")
})

test_that("new_gdpar_coef rejects non-integer or non-positive p", {
  theta_df <- make_theta_df(p = 1L)
  expect_error(new_gdpar_coef(theta_ref = theta_df, p = 0L),
               class = "gdpar_internal_error")
  expect_error(new_gdpar_coef(theta_ref = theta_df, p = -1L),
               class = "gdpar_internal_error")
  expect_error(new_gdpar_coef(theta_ref = theta_df, p = 1.5),
               class = "gdpar_internal_error")
  expect_error(new_gdpar_coef(theta_ref = theta_df, p = c(1L, 2L)),
               class = "gdpar_internal_error")
})

test_that("new_gdpar_coef rejects theta_ref with wrong cols", {
  bad <- data.frame(k = 1L, mean = 0.1)
  expect_error(new_gdpar_coef(theta_ref = bad, p = 1L),
               class = "gdpar_internal_error")
})

test_that("new_gdpar_coef rejects theta_ref with wrong nrow", {
  theta_df <- make_theta_df(p = 2L)
  expect_error(new_gdpar_coef(theta_ref = theta_df, p = 3L),
               class = "gdpar_internal_error")
})

test_that("new_gdpar_coef rejects theta_ref with non-canonical k", {
  bad <- data.frame(k = c(1L, 3L), mean = c(0, 0),
                    q05 = c(0, 0), q50 = c(0, 0), q95 = c(0, 0))
  expect_error(new_gdpar_coef(theta_ref = bad, p = 2L),
               class = "gdpar_internal_error")
})

test_that("new_gdpar_coef rejects slot of wrong length", {
  theta_df <- make_theta_df(p = 2L)
  a_df <- make_term_df()
  expect_error(
    new_gdpar_coef(theta_ref = theta_df, a = list(a_df), p = 2L),
    class = "gdpar_internal_error"
  )
})

test_that("new_gdpar_coef rejects non-data.frame slot entry", {
  theta_df <- make_theta_df(p = 2L)
  expect_error(
    new_gdpar_coef(theta_ref = theta_df,
                   a = list("notadf", NULL), p = 2L),
    class = "gdpar_internal_error"
  )
})

test_that("new_gdpar_coef rejects slot entry missing required cols", {
  theta_df <- make_theta_df(p = 1L)
  bad_df <- data.frame(term = "x1", mean = 0)
  expect_error(
    new_gdpar_coef(theta_ref = theta_df,
                   a = list(bad_df), p = 1L),
    class = "gdpar_internal_error"
  )
})

# ---- S3 methods ----------------------------------------------------

test_that("print.gdpar_coef returns invisibly and supports three levels", {
  theta_df <- make_theta_df(p = 2L)
  a_df <- make_term_df()
  obj <- new_gdpar_coef(theta_ref = theta_df,
                        a = list(a_df, NULL), p = 2L)
  expect_invisible(print(obj))
  out_g <- capture.output(print(obj, level = "global"))
  out_c <- capture.output(print(obj, level = "coord"))
  out_f <- capture.output(print(obj, level = "full"))
  expect_true(any(grepl("<gdpar_coef>", out_g)))
  expect_true(any(grepl("coord k = 1", out_c)))
  expect_true(any(grepl("coord k = 1", out_f)))
  expect_false(any(grepl("coord k = 1", out_g)))
})

test_that("print.gdpar_coef rejects invalid level", {
  theta_df <- make_theta_df(p = 1L)
  obj <- new_gdpar_coef(theta_ref = theta_df, p = 1L)
  expect_error(print(obj, level = "nope"))
})

test_that("summary.gdpar_coef counts active components correctly", {
  theta_df <- make_theta_df(p = 3L)
  a_df <- make_term_df()
  W_df <- make_W_df()
  obj <- new_gdpar_coef(
    theta_ref = theta_df,
    a = list(a_df, NULL, a_df),
    b = NULL,
    W = list(NULL, W_df, NULL),
    p = 3L
  )
  s <- summary(obj)
  expect_equal(s$p, 3L)
  expect_equal(unname(s$n_active["a"]), 2L)
  expect_equal(unname(s$n_active["b"]), 0L)
  expect_equal(unname(s$n_active["W"]), 1L)
})

test_that("format.gdpar_coef returns a one-line summary", {
  theta_df <- make_theta_df(p = 2L)
  obj <- new_gdpar_coef(theta_ref = theta_df, p = 2L)
  s <- format(obj)
  expect_type(s, "character")
  expect_length(s, 1L)
  expect_match(s, "^<gdpar_coef>")
})

test_that("as.data.frame.gdpar_coef returns the expected long-tidy schema", {
  theta_df <- make_theta_df(p = 2L)
  a_df <- make_term_df(n_terms = 2L)
  W_df <- make_W_df(basis_dim = 2L, d_x = 2L)
  obj <- new_gdpar_coef(
    theta_ref = theta_df,
    a = list(a_df, NULL),
    W = list(W_df, NULL),
    p = 2L
  )
  df <- as.data.frame(obj)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("component", "g", "k", "identifier", "x_name",
                     "mean", "q05", "q50", "q95"))
  expect_true(all(c("theta_ref", "a", "W") %in% df$component))
  expect_equal(sum(df$component == "theta_ref"), 2L)
  expect_equal(sum(df$component == "a"), 2L)
  expect_equal(sum(df$component == "W"), 4L)
  expect_true(all(is.na(df$x_name[df$component == "theta_ref"])))
  expect_true(all(!is.na(df$x_name[df$component == "W"])))
  expect_true(all(is.na(df$g)))
})

test_that("as.data.frame.gdpar_coef emits g and hyperparam rows under grouping", {
  arr <- array(rnorm(40L * 3L * 2L), dim = c(40L, 3L, 2L))
  theta_grouped <- build_coef_theta_ref_df_grouped(arr, J_groups = 3L,
                                                    p = 2L)
  mu_df <- data.frame(k = 1:2, mean = c(0.1, 0.2),
                      q05 = c(-0.1, 0), q50 = c(0.1, 0.2),
                      q95 = c(0.3, 0.4))
  sigma_df <- data.frame(k = 1:2, mean = c(1.0, 1.1),
                         q05 = c(0.5, 0.6), q50 = c(1.0, 1.1),
                         q95 = c(1.5, 1.6))
  obj <- new_gdpar_coef(
    theta_ref       = theta_grouped,
    a               = NULL,
    b               = NULL,
    W               = NULL,
    p               = 2L,
    mu_theta_ref    = mu_df,
    sigma_theta_ref = sigma_df,
    J_groups        = 3L,
    group_levels    = c("a", "b", "c")
  )
  df <- as.data.frame(obj)
  expect_named(df, c("component", "g", "k", "identifier", "x_name",
                     "mean", "q05", "q50", "q95"))
  expect_equal(sum(df$component == "theta_ref"), 6L)
  expect_equal(sum(df$component == "mu_theta_ref"), 2L)
  expect_equal(sum(df$component == "sigma_theta_ref"), 2L)
  expect_true(all(!is.na(df$g[df$component == "theta_ref"])))
  expect_true(all(is.na(df$g[df$component %in%
                                c("mu_theta_ref", "sigma_theta_ref")])))
})

test_that("as.data.frame.gdpar_coef handles empty (no components) gracefully", {
  theta_df <- make_theta_df(p = 1L)
  obj <- new_gdpar_coef(theta_ref = theta_df, p = 1L)
  df <- as.data.frame(obj)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_equal(df$component, "theta_ref")
})

test_that("as.data.frame.gdpar_coef preserves numeric values bit-exact for theta_ref", {
  theta_df <- make_theta_df(p = 2L)
  obj <- new_gdpar_coef(theta_ref = theta_df, p = 2L)
  df <- as.data.frame(obj)
  expect_equal(df$mean[df$component == "theta_ref"], theta_df$mean)
  expect_equal(df$q05[df$component == "theta_ref"], theta_df$q05)
})

# ---- count_active_coords helper ------------------------------------

test_that("count_active_coords handles NULL and empty correctly", {
  expect_equal(count_active_coords(NULL), 0L)
  expect_equal(count_active_coords(list(NULL, NULL)), 0L)
  df <- data.frame(term = "x", mean = 1, q05 = 0, q50 = 1, q95 = 2)
  expect_equal(count_active_coords(list(df, NULL, df)), 2L)
})

test_that("count_active_coords excludes data.frames with zero rows", {
  empty <- data.frame(term = character(0), mean = numeric(0),
                      q05 = numeric(0), q50 = numeric(0),
                      q95 = numeric(0))
  expect_equal(count_active_coords(list(empty, NULL)), 0L)
})

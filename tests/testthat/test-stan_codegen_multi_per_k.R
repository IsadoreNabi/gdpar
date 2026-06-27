# Tests for Phase H.2: per-k CP/NCP wiring in the multivariate
# Stan template. Coverage: codegen branches (uniform NCP, uniform
# CP, mixed), bit-exact paridad with H.1 NCP-only template, Stan
# syntax check for each branch, codegen validation errors.

prior_min <- function() gdpar_prior()


# ---- bit-exact paridad: uniform NCP equals H.1 codegen ----------

test_that("codegen multi uniform NCP via scalar cp_a is byte-identical to cp_a_per_k all FALSE", {
  s1 <- gdpar:::generate_stan_code_multi(prior_min(),
                                         cp_a = FALSE, cp_W = FALSE)
  s2 <- gdpar:::generate_stan_code_multi(prior_min(),
                                         cp_a_per_k = c(FALSE, FALSE))
  s3 <- gdpar:::generate_stan_code_multi(prior_min(),
                                         cp_a_per_k = c(FALSE, FALSE, FALSE))
  expect_identical(s1, s2)
  expect_identical(s2, s3)
})

test_that("codegen multi uniform NCP retains the H.1 sampling statement", {
  s <- gdpar:::generate_stan_code_multi(prior_min(),
                                        cp_a = FALSE, cp_W = FALSE)
  expect_match(s, "a_raw ~ normal\\(0, 1\\);")
  expect_false(grepl("segment\\(a_raw", s))
  expect_false(grepl("cp_a_per_k_data", s))
  expect_match(s, "a_coef\\[k\\]\\[j\\] = a_raw\\[a_raw_offset\\[k\\] \\+ j\\] \\* sigma_a\\[k\\];")
})


# ---- uniform CP branch ------------------------------------------

test_that("codegen multi uniform CP wires segment-based per-k priors", {
  s <- gdpar:::generate_stan_code_multi(prior_min(),
                                        cp_a = TRUE, cp_W = FALSE)
  expect_match(s, "segment\\(a_raw, a_raw_offset\\[k\\] \\+ 1, J_a_free\\[k\\]\\)")
  expect_match(s, "~ normal\\(0, sigma_a\\[k\\]\\);")
  expect_false(grepl("a_raw ~ normal", s))
  expect_false(grepl("cp_a_per_k_data", s))
  # In CP, a_coef should NOT carry the * sigma_a[k] scale in tp.
  expect_false(grepl("a_coef\\[k\\]\\[j\\] = a_raw\\[a_raw_offset\\[k\\] \\+ j\\] \\* sigma_a",
                     s))
  expect_match(s, "a_coef\\[k\\]\\[j\\] = a_raw\\[a_raw_offset\\[k\\] \\+ j\\];")
})

test_that("codegen multi uniform CP via cp_a_per_k all TRUE matches scalar cp_a = TRUE", {
  s_scalar <- gdpar:::generate_stan_code_multi(prior_min(),
                                               cp_a = TRUE, cp_W = FALSE)
  s_vec    <- gdpar:::generate_stan_code_multi(prior_min(),
                                               cp_a_per_k = c(TRUE, TRUE))
  expect_identical(s_scalar, s_vec)
})


# ---- mixed branch ------------------------------------------------

test_that("codegen multi mixed cp_a_per_k declares cp_a_per_k_data and ternary", {
  s <- gdpar:::generate_stan_code_multi(prior_min(),
                                        cp_a_per_k = c(TRUE, FALSE))
  expect_match(s, "array\\[p\\] int<lower=0, upper=1> cp_a_per_k_data;")
  expect_match(s, "real scale_k = \\(cp_a_per_k_data\\[k\\] == 1\\) \\? 1.0 : sigma_a\\[k\\];")
  expect_match(s, "if \\(cp_a_per_k_data\\[k\\] == 1\\)")
  expect_match(s, "~ normal\\(0, sigma_a\\[k\\]\\);")
  expect_match(s, "~ normal\\(0, 1\\);")
  expect_false(grepl("a_raw ~ normal\\(0, 1\\);", s))
})

test_that("codegen multi mixed branch declarations are absent in uniform branches", {
  s_uniform_ncp <- gdpar:::generate_stan_code_multi(prior_min(),
                                                    cp_a_per_k = c(FALSE, FALSE))
  s_uniform_cp <- gdpar:::generate_stan_code_multi(prior_min(),
                                                   cp_a_per_k = c(TRUE, TRUE))
  for (s in list(s_uniform_ncp, s_uniform_cp)) {
    expect_false(grepl("cp_a_per_k_data", s))
    expect_false(grepl("scale_k", s))
  }
})


# ---- validation --------------------------------------------------

test_that("generate_stan_code_multi rejects cp_a_per_k with NA", {
  expect_error(
    gdpar:::generate_stan_code_multi(prior_min(),
                                     cp_a_per_k = c(TRUE, NA)),
    class = "gdpar_input_error"
  )
})

test_that("generate_stan_code_multi rejects non-logical cp_a_per_k", {
  expect_error(
    gdpar:::generate_stan_code_multi(prior_min(),
                                     cp_a_per_k = c(1, 0)),
    class = "gdpar_input_error"
  )
})

test_that("generate_stan_code_multi rejects empty cp_a_per_k", {
  expect_error(
    gdpar:::generate_stan_code_multi(prior_min(),
                                     cp_a_per_k = logical(0)),
    class = "gdpar_input_error"
  )
})


# ---- helper unit tests ------------------------------------------

test_that("generate_a_blocks_multi: all FALSE returns empty data_decl", {
  bl <- gdpar:::generate_a_blocks_multi(c(FALSE, FALSE))
  expect_identical(bl$data_decl, "")
  expect_match(bl$model_block, "a_raw ~ normal\\(0, 1\\);")
})

test_that("generate_a_blocks_multi: all TRUE returns empty data_decl with CP block", {
  bl <- gdpar:::generate_a_blocks_multi(c(TRUE, TRUE, TRUE))
  expect_identical(bl$data_decl, "")
  expect_match(bl$model_block, "segment\\(a_raw")
  expect_match(bl$model_block, "normal\\(0, sigma_a\\[k\\]\\)")
})

test_that("generate_a_blocks_multi: mixed returns array decl + ternary", {
  bl <- gdpar:::generate_a_blocks_multi(c(TRUE, FALSE))
  expect_match(bl$data_decl,
               "\\n  array\\[p\\] int<lower=0, upper=1> cp_a_per_k_data;")
  expect_match(bl$tp_block, "scale_k = \\(cp_a_per_k_data\\[k\\] == 1\\)")
})


# ---- env-gated: Stan syntax check on each branch ----------------

test_that("Stan template H.2 compiles for uniform NCP, uniform CP, and mixed", {
  skip_if_not_installed("cmdstanr")
  for (pat in list(c(FALSE, FALSE),
                   c(TRUE,  TRUE),
                   c(TRUE,  FALSE),
                   c(FALSE, TRUE))) {
    src <- gdpar:::generate_stan_code_multi(prior_min(),
                                            cp_a_per_k = pat,
                                            cp_W = FALSE)
    tmpfile <- tempfile(fileext = ".stan")
    writeLines(src, tmpfile)
    mod <- cmdstanr::cmdstan_model(tmpfile, compile = FALSE)
    expect_true(mod$check_syntax(pedantic = FALSE, quiet = TRUE),
                info = sprintf("pattern: %s",
                               paste(pat, collapse = ", ")))
  }
})


# ---- env-gated: end-to-end mixed per-k via gdpar() --------------

test_that("gdpar() multi with parametrization_aggregation='per_k' wires cp_a_per_k_data", {
  testthat::skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
    "Set GDPAR_F6_FIT=1 to run the multi per-k end-to-end smoke (compiles a Stan model; takes minutes).")
  skip_if_not_installed("cmdstanr")
  set.seed(42)
  n <- 80L
  x1 <- rnorm(n); x2 <- rnorm(n)
  x1_c <- x1 - mean(x1); x2_c <- x2 - mean(x2)
  # Coord 1: high-info (CP expected). Coord 2: low-info (NCP expected).
  y1 <- 0.5 + 1.0 * x1_c - 0.8 * x2_c + rnorm(n, sd = 0.1)
  y2 <- 0.4 + 0.02 * x1_c - 0.01 * x2_c + rnorm(n, sd = 2.0)
  df <- data.frame(x1 = x1, x2 = x2)
  df$y <- cbind(y1, y2)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  fit <- suppressMessages(suppressWarnings(
    gdpar(
      formula = y ~ x1 + x2, family = fam, amm = spec, data = df,
      chains = 2L, iter_warmup = 200L, iter_sampling = 200L,
      refresh = 0L, verbose = FALSE, seed = 42L,
      parametrization = "auto",
      parametrization_aggregation = "per_k",
      skip_id_check = TRUE
    )
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_identical(fit$parametrization$meta$aggregation, "per_k")
  expect_s3_class(fit$parametrization$report,
                  "gdpar_preflight_report")
  # cp_a per_k decisions are vectors; cp_a global is NA under per_k
  # because there is no scalar interpretation.
  expect_true(is.na(fit$parametrization$cp_a))
})

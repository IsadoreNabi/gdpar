# ============================================================
# gdpar / inst/benchmarks/quick_smoke.R
# Fast end-to-end smoke test of the public API.
# Designed to finish in ~2-5 minutes on a modern multi-core CPU.
# Exercises every public function on small synthetic datasets and
# reports pass/fail per check. Intended as a first sanity check
# before running the heavier validation scripts.
# ============================================================

if (!exists("bench_section", mode = "function")) {
  helpers_path <- local({
    candidates <- character(0)
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      ctx <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                      error = function(e) "")
      if (nzchar(ctx) && file.exists(ctx)) {
        candidates <- c(candidates,
                        file.path(dirname(ctx), "helpers.R"))
      }
    }
    sourced <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(sourced) && nzchar(sourced)) {
      candidates <- c(candidates,
                      file.path(dirname(normalizePath(sourced,
                                                     mustWork = FALSE)),
                                "helpers.R"))
    }
    candidates <- c(candidates,
                    file.path(getwd(), "helpers.R"),
                    file.path(getwd(), "inst", "benchmarks",
                              "helpers.R"))
    inst_path <- system.file("benchmarks", "helpers.R",
                             package = "gdpar")
    if (nzchar(inst_path)) candidates <- c(candidates, inst_path)
    found <- ""
    for (p in candidates) {
      if (nzchar(p) && file.exists(p)) { found <- p; break }
    }
    if (!nzchar(found)) {
      stop("Could not locate helpers.R. Set the working directory ",
           "to the gdpar package root or to its inst/benchmarks ",
           "subdirectory, or install gdpar so that system.file() ",
           "resolves it.")
    }
    found
  })
  source(helpers_path)
}

bench_section("gdpar quick smoke test")

if (!bench_require("gdpar")) stop("gdpar must be installed first.")
library(gdpar)
have_cmdstan <- bench_have_cmdstan()

n_checks <- 0L
n_passed <- 0L
record <- function(ok, msg) {
  n_checks <<- n_checks + 1L
  if (isTRUE(ok)) {
    n_passed <<- n_passed + 1L
    bench_pass(msg)
  } else {
    bench_fail(msg)
  }
}

# ----- A. Constructors -----

bench_subsection("A. Constructors")

ok <- tryCatch({
  fam <- gdpar_family("gaussian")
  inherits(fam, "gdpar_family") && fam$name == "gaussian"
}, error = function(e) FALSE)
record(ok, "gdpar_family('gaussian') constructs a gdpar_family object")

ok <- tryCatch({
  fam <- gdpar_family("poisson")
  inherits(fam, "gdpar_family") && fam$name == "poisson"
}, error = function(e) FALSE)
record(ok, "gdpar_family('poisson') constructs a gdpar_family object")

ok <- tryCatch({
  fam <- gdpar_family("neg_binomial_2")
  inherits(fam, "gdpar_family") && fam$has_dispersion
}, error = function(e) FALSE)
record(ok, "gdpar_family('neg_binomial_2') has dispersion = TRUE")

ok <- tryCatch({
  fam <- gdpar_family("bernoulli")
  inherits(fam, "gdpar_family") && fam$link == "logit"
}, error = function(e) FALSE)
record(ok, "gdpar_family('bernoulli') uses logit link by default")

ok <- tryCatch({
  inherits(gdpar_family_custom(
    name = "qsmoke_custom", link = "log", did_holds = TRUE,
    did_condition = NA_character_,
    stan_loglik_block = "target += normal_lpdf(y[i] | eta[i], 1);",
    did_reference = "smoke test"
  ), "gdpar_family")
}, error = function(e) FALSE)
record(ok, "gdpar_family_custom constructs a gdpar_family object")

ok <- tryCatch({
  wb <- W_basis(type = "polynomial", degree = 2)
  inherits(wb, "W_basis") && length(wb$evaluator(0.5)) == 2L
}, error = function(e) FALSE)
record(ok, "W_basis polynomial degree 2 evaluator returns length 2 at theta = 0.5")

ok <- tryCatch({
  spec <- amm_spec(a = ~ x1 + x2, b = ~ x1, W = W_basis("polynomial", 1))
  inherits(spec, "amm_spec") && spec$level == 2L
}, error = function(e) FALSE)
record(ok, "amm_spec with a, b, W active records Level 2")

ok <- tryCatch({
  pr <- gdpar_prior()
  inherits(pr, "gdpar_prior")
}, error = function(e) FALSE)
record(ok, "gdpar_prior() returns a gdpar_prior object with defaults")

# ----- B. Identifiability diagnostic -----

bench_subsection("B. Identifiability diagnostic")

ok <- tryCatch({
  df <- make_gauss_clean(n = 200, seed = 42)
  rep <- gdpar_check_identifiability(amm_spec(a = ~ x1 + x2), df,
                                      formula_rhs = ~ x1 + x2)
  isTRUE(rep$passed)
}, error = function(e) FALSE)
record(ok, "Independent additive basis passes Gram diagnostic")

ok <- tryCatch({
  df <- make_gauss_clean(n = 200, seed = 42)
  df$x_dup <- df$x1
  rep <- gdpar_check_identifiability(
    amm_spec(a = ~ x1 + x_dup), df, formula_rhs = ~ x1 + x_dup
  )
  isFALSE(rep$passed)
}, error = function(e) FALSE)
record(ok, "Duplicated covariate fails Gram diagnostic")

ok <- tryCatch({
  rep <- gdpar_check_identifiability(amm_spec(), data.frame(x = rnorm(50)))
  isTRUE(rep$passed)
}, error = function(e) FALSE)
record(ok, "Level 0 identifiability check passes trivially")

# ----- C. Stan code generation (does not require cmdstan) -----

bench_subsection("C. Stan code generation")

ok <- tryCatch({
  src <- gdpar:::generate_stan_code(gdpar_prior())
  is.character(src) && nchar(src) > 100L &&
    !grepl("{{", src, fixed = TRUE)
}, error = function(e) FALSE)
record(ok, "generate_stan_code substitutes all placeholders")

ok <- tryCatch({
  src <- gdpar:::generate_stan_code(
    gdpar_prior(theta_ref = "normal(2, 5)")
  )
  grepl("normal(2, 5)", src, fixed = TRUE)
}, error = function(e) FALSE)
record(ok, "generate_stan_code uses user-supplied prior strings")

ok <- tryCatch({
  src <- gdpar:::generate_stan_code(gdpar_prior())
  path <- gdpar:::write_stan_to_tempfile(src)
  startsWith(path, tempdir()) && file.exists(path)
}, error = function(e) FALSE)
record(ok, "write_stan_to_tempfile writes inside tempdir")

# ----- D. End-to-end fit (Gaussian) -----

bench_subsection("D. End-to-end fit (Gaussian, AMM Level 1)")

if (have_cmdstan) {
  fit_gauss <- tryCatch({
    df <- make_gauss_clean(n = 200, seed = 42)
    gdpar(
      formula       = y ~ x1 + x2,
      family        = gdpar_family("gaussian"),
      amm           = amm_spec(a = ~ x1 + x2),
      data          = df,
      iter_warmup   = 300,
      iter_sampling = 300,
      chains        = 2,
      refresh       = 0,
      verbose       = FALSE
    )
  }, error = function(e) {
    record(FALSE, sprintf("Gaussian fit raised an error: %s",
                          conditionMessage(e)))
    NULL
  })
  if (!is.null(fit_gauss)) {
    record(inherits(fit_gauss, "gdpar_fit"),
           "gdpar(gaussian, AMM Level 1) returns a gdpar_fit")
    record(isTRUE(fit_gauss$diagnostics$converged) ||
             fit_gauss$diagnostics$rhat_max <= 1.10,
           sprintf("Gaussian fit converged (rhat_max = %.3f)",
                   fit_gauss$diagnostics$rhat_max))
    record(!is.null(coef(fit_gauss)$theta_ref),
           "coef(fit) returns theta_ref")
    pr <- predict(fit_gauss, summary = "mean_se")
    record(is.data.frame(pr) && nrow(pr) == 200L,
           "predict(summary='mean_se') returns one row per observation")
    pr_resp <- predict(fit_gauss, type = "response", summary = "mean_se")
    record(is.data.frame(pr_resp) && nrow(pr_resp) == 200L,
           "predict(type='response') returns one row per observation")
    s <- summary(fit_gauss)
    record(is.data.frame(s) && nrow(s) >= 1L,
           "summary(fit) returns a data.frame")
    print_ok <- tryCatch({
      invisible(capture.output(print(fit_gauss)))
      TRUE
    }, error = function(e) FALSE)
    record(print_ok, "print(fit) executes without error")
  }
} else {
  bench_skip("Gaussian end-to-end fit skipped (cmdstan unavailable)")
}

# ----- E. End-to-end fit (Poisson) -----

bench_subsection("E. End-to-end fit (Poisson, AMM Level 1)")

if (have_cmdstan) {
  fit_pois <- tryCatch({
    df <- make_poisson_clean(n = 200, seed = 42)
    gdpar(
      formula       = y ~ x1 + x2,
      family        = gdpar_family("poisson"),
      amm           = amm_spec(a = ~ x1 + x2),
      data          = df,
      iter_warmup   = 300,
      iter_sampling = 300,
      chains        = 2,
      refresh       = 0,
      verbose       = FALSE
    )
  }, error = function(e) {
    record(FALSE, sprintf("Poisson fit raised an error: %s",
                          conditionMessage(e)))
    NULL
  })
  if (!is.null(fit_pois)) {
    record(inherits(fit_pois, "gdpar_fit"),
           "gdpar(poisson, AMM Level 1) returns a gdpar_fit")
    record(isTRUE(fit_pois$diagnostics$converged) ||
             fit_pois$diagnostics$rhat_max <= 1.10,
           sprintf("Poisson fit converged (rhat_max = %.3f)",
                   fit_pois$diagnostics$rhat_max))
  }
} else {
  bench_skip("Poisson end-to-end fit skipped (cmdstan unavailable)")
}

# ----- F. AMM Level 2 fit -----

bench_subsection("F. End-to-end fit (Gaussian AMM Level 2 with active b)")

if (have_cmdstan) {
  fit_l2 <- tryCatch({
    df <- make_gauss_amm_level2(n = 300, seed = 13)
    gdpar(
      formula       = y ~ x1 + x2 + x3 + x4,
      family        = gdpar_family("gaussian"),
      amm           = amm_spec(a = ~ x1 + x2, b = ~ x3 + x4),
      data          = df,
      iter_warmup   = 400,
      iter_sampling = 400,
      chains        = 2,
      refresh       = 0,
      verbose       = FALSE,
      anchor        = "empirical_y"
    )
  }, error = function(e) {
    record(FALSE, sprintf("AMM Level 2 fit raised an error: %s",
                          conditionMessage(e)))
    NULL
  })
  if (!is.null(fit_l2)) {
    record(inherits(fit_l2, "gdpar_fit"),
           "gdpar(gaussian, AMM Level 2 with active b) returns a gdpar_fit")
    record(fit_l2$diagnostics$rhat_max <= 1.10,
           sprintf("AMM Level 2 fit converged (rhat_max = %.3f)",
                   fit_l2$diagnostics$rhat_max))
  }
} else {
  bench_skip("AMM Level 2 fit skipped (cmdstan unavailable)")
}

# ----- G. Final summary -----

bench_section("Quick smoke summary")
cat(sprintf("Checks passed: %d / %d\n", n_passed, n_checks))
if (n_passed < n_checks) {
  cat("Some checks failed. Inspect the [FAIL] lines above.\n")
} else {
  cat("All quick smoke checks passed.\n")
}

invisible(NULL)

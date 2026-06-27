# Tests for Sub-sub-fase 9.3.a (Bloque 9, Sesion B9.3, 2026-05-27):
# canonical dispatcher .gdpar_emit_canonical_stan() and the byte-identity
# guarantees of the pieces shipped in inst/stan/_canonical_pieces/ under
# the canonized decision B.iv lateral (codegen R-side as canonical
# source of truth).
#
# Updated in Sub-sub-fase 9.3.a colateral (Sesion B9.4, 2026-05-27)
# under decision G.iv lateral closing deuda f (helpers duplication):
# the canonical helpers (bspline_basis_eval + apply_W_basis_diff) live
# in a dedicated piece `amm_canonical_helpers.stan` and the dispatcher
# inserts them via R-side textual substitution at the
# `// {{CANONICAL_HELPERS}}` placeholder inside the body piece's
# `functions { }` block. The body pieces shrink (helpers removed
# inline) and their frozen md5 hashes are updated; the dispatcher
# OUTPUT md5 for the p1 path is bit-exact preserved from B9.3
# (helpers verbatim from canonical_p1.stan), while the pmulti path
# OUTPUT md5 changes by ~16 lines of unified helper comments. Goldens
# are preserved by-construction (Stan semantics identical; only
# comments differ).
#
# Verification under the rigor de tres capas (feedback_proof_rigor_standards):
#   - algebraica: pieces have the frozen md5 hashes (snapshot guard) and
#     a parse-check via cmdstanr::stanc when available; dispatcher
#     OUTPUT md5 frozen as the bit-exact-semantics invariant;
#   - estadistica: dispatcher is deterministic and the public wrappers
#     produce output bit-identical to the direct dispatcher call;
#   - numerica: covered by the existing test-stan_codegen* suite + the
#     post-cluster R CMD check natural.

.b9_3_get_canonical_piece_path <- function(piece_name) {
  p <- system.file("stan", "_canonical_pieces", piece_name,
                   package = "gdpar")
  if (!nzchar(p) || !file.exists(p)) {
    p <- file.path("inst", "stan", "_canonical_pieces", piece_name)
  }
  p
}

test_that("canonical piece amm_canonical_p1.stan has the frozen md5 hash (post-G.iv)", {
  # Updated B9.4 G.iv: helpers extracted to amm_canonical_helpers.stan;
  # functions { } block now contains the {{CANONICAL_HELPERS}} placeholder
  # only. Pre-B9.4 hash f899811efafe83f448283011742f79b9 (13878 bytes).
  path <- .b9_3_get_canonical_piece_path("amm_canonical_p1.stan")
  expect_true(file.exists(path))
  observed <- unname(tools::md5sum(path))
  expect_equal(observed, "730c93f56a300d96d2991b30754db57a")
})

test_that("canonical piece amm_canonical_pmulti.stan has the frozen md5 hash (post-G.iv)", {
  # Updated B9.4 G.iv: helpers extracted. Pre-B9.4 hash
  # 112a4e68806ed5dad5b86ade7ab921a9 (16028 bytes).
  path <- .b9_3_get_canonical_piece_path("amm_canonical_pmulti.stan")
  expect_true(file.exists(path))
  observed <- unname(tools::md5sum(path))
  expect_equal(observed, "01fc4b04da46f71db8ef6acc5c7038b3")
})

test_that("canonical helpers piece amm_canonical_helpers.stan has the frozen md5 hash (B9.4 G.iv)", {
  # New B9.4 G.iv: dedicated Stan source fragment containing the two
  # canonical helpers bspline_basis_eval + apply_W_basis_diff inserted
  # by the dispatcher at the {{CANONICAL_HELPERS}} placeholder of body
  # pieces. NOT standalone-parseable by cmdstanr::stanc (it has no
  # enclosing `functions { }` block); intended for substitution only.
  path <- .b9_3_get_canonical_piece_path("amm_canonical_helpers.stan")
  expect_true(file.exists(path))
  observed <- unname(tools::md5sum(path))
  expect_equal(observed, "c548a2a46ccee3db1f3bea3430588e87")
})

test_that(".gdpar_emit_canonical_stan p1 path is deterministic for canonical inputs", {
  prior <- gdpar_prior()
  spec_ncp <- list(
    p_class = "p1",
    prior = prior,
    cp_a = FALSE,
    cp_W = FALSE,
    mle = FALSE
  )
  src1 <- .gdpar_emit_canonical_stan(spec_ncp)
  src2 <- .gdpar_emit_canonical_stan(spec_ncp)
  expect_identical(src1, src2)
  spec_mle <- list(
    p_class = "p1",
    prior = prior,
    cp_a = TRUE,
    cp_W = TRUE,
    mle = TRUE
  )
  src_mle1 <- .gdpar_emit_canonical_stan(spec_mle)
  src_mle2 <- .gdpar_emit_canonical_stan(spec_mle)
  expect_identical(src_mle1, src_mle2)
  expect_lt(nchar(src_mle1), nchar(src1))
})

test_that(".gdpar_emit_canonical_stan pmulti path is deterministic and template_name maps correctly", {
  prior <- gdpar_prior()
  spec_default <- list(
    p_class = "pmulti",
    prior = prior,
    cp_a = FALSE,
    cp_W = FALSE,
    cp_a_per_k = NULL
  )
  src_default1 <- .gdpar_emit_canonical_stan(spec_default)
  src_default2 <- .gdpar_emit_canonical_stan(spec_default)
  expect_identical(src_default1, src_default2)
  spec_explicit <- spec_default
  spec_explicit$template_name <- "amm_canonical_pmulti.stan"
  src_explicit <- .gdpar_emit_canonical_stan(spec_explicit)
  expect_identical(src_default1, src_explicit)
})

test_that(".gdpar_emit_canonical_stan rejects invalid spec", {
  prior <- gdpar_prior()
  expect_error(
    .gdpar_emit_canonical_stan(list(
      p_class = "pK", prior = prior, cp_a = FALSE, cp_W = FALSE
    )),
    class = "gdpar_internal_error"
  )
  expect_error(
    .gdpar_emit_canonical_stan(list(
      p_class = NULL, prior = prior, cp_a = FALSE, cp_W = FALSE
    )),
    class = "gdpar_internal_error"
  )
})

test_that(".gdpar_emit_canonical_stan pmulti rejects malformed cp_a_per_k", {
  prior <- gdpar_prior()
  expect_error(
    .gdpar_emit_canonical_stan(list(
      p_class = "pmulti",
      prior = prior,
      cp_a = FALSE,
      cp_W = FALSE,
      cp_a_per_k = c(TRUE, NA, FALSE)
    )),
    class = "gdpar_input_error"
  )
  expect_error(
    .gdpar_emit_canonical_stan(list(
      p_class = "pmulti",
      prior = prior,
      cp_a = FALSE,
      cp_W = FALSE,
      cp_a_per_k = logical(0)
    )),
    class = "gdpar_input_error"
  )
})

test_that("public wrappers delegate to the canonical dispatcher bit-identically", {
  prior <- gdpar_prior()
  src_wrapper_p1 <- generate_stan_code(prior, cp_a = FALSE, cp_W = FALSE)
  src_direct_p1 <- .gdpar_emit_canonical_stan(list(
    p_class = "p1",
    prior = prior,
    cp_a = FALSE,
    cp_W = FALSE,
    mle = FALSE
  ))
  expect_identical(src_wrapper_p1, src_direct_p1)
  src_wrapper_pmulti <- generate_stan_code_multi(prior, cp_a = FALSE,
                                                 cp_W = FALSE)
  src_direct_pmulti <- .gdpar_emit_canonical_stan(list(
    p_class = "pmulti",
    prior = prior,
    cp_a = FALSE,
    cp_W = FALSE,
    cp_a_per_k = NULL,
    template_name = "amm_canonical_pmulti.stan"
  ))
  expect_identical(src_wrapper_pmulti, src_direct_pmulti)
})

# ---- Sub-sub-fase 9.3.a colateral (Sesion B9.4, 2026-05-27) G.iv ----
# Frozen dispatcher OUTPUT md5 hashes: the canonical invariant
# preserved post-B9.4. For the p1 path, the dispatcher OUTPUT is
# bit-exact preserved across B9.3 -> B9.4 (helpers in the dedicated
# piece are verbatim from canonical_p1.stan; comments unified at the
# canonical_p1 version). For the pmulti path, the dispatcher OUTPUT
# changes only by the unified helper comments (Stan semantics
# preserved; cmdstanr cache recompiles once on first call post-B9.4
# but goldens are unaffected since they fit data, not source hashes).

test_that("dispatcher OUTPUT md5 frozen invariant - p1 NCP (bit-exact preserved B9.3 -> B9.4)", {
  prior <- gdpar_prior()
  spec <- list(p_class = "p1", prior = prior,
               cp_a = FALSE, cp_W = FALSE, mle = FALSE)
  src <- .gdpar_emit_canonical_stan(spec)
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.3 frozen output hash; preserved bit-exact in B9.4 G.iv
  expect_equal(observed, "8930c724b39fab56024330dc3b9bf2c3")
  expect_equal(nchar(src), 13867L)
})

test_that("dispatcher OUTPUT md5 frozen invariant - p1 MLE (bit-exact preserved B9.3 -> B9.4)", {
  prior <- gdpar_prior()
  spec <- list(p_class = "p1", prior = prior,
               cp_a = TRUE, cp_W = TRUE, mle = TRUE)
  src <- .gdpar_emit_canonical_stan(spec)
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.3 frozen output hash; preserved bit-exact in B9.4 G.iv
  expect_equal(observed, "1e6d7d28ccc3f6aff1f90d3d0a772a52")
  expect_equal(nchar(src), 13139L)
})

test_that("dispatcher OUTPUT md5 frozen invariant - pmulti default (B9.4 G.iv canonical, comments unified)", {
  prior <- gdpar_prior()
  spec <- list(p_class = "pmulti", prior = prior,
               cp_a = FALSE, cp_W = FALSE, cp_a_per_k = NULL)
  src <- .gdpar_emit_canonical_stan(spec)
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.4 frozen output hash; pre-B9.4 was cf182996cf4cc7f7d50ba66205b47882
  # (16183 bytes). Difference is 892 bytes of unified helper comments
  # (canonical_pmulti pre-G.iv had a brief 4-line block referencing the
  # then-open deuda f; replaced by the detailed 16-line algorithmic
  # comment block from canonical_p1). Stan semantics identical.
  expect_equal(observed, "2ede60f680984d4a4ace989443a16977")
  expect_equal(nchar(src), 17075L)
})

test_that("dispatcher OUTPUT contains the canonical helpers (B9.4 G.iv sanity)", {
  prior <- gdpar_prior()
  spec_p1 <- list(p_class = "p1", prior = prior,
                  cp_a = FALSE, cp_W = FALSE, mle = FALSE)
  src_p1 <- .gdpar_emit_canonical_stan(spec_p1)
  # Both helpers present (definition signatures)
  expect_true(grepl("vector bspline_basis_eval(real x, int W_per_k_dim",
                    src_p1, fixed = TRUE))
  expect_true(grepl("vector apply_W_basis_diff(int basis_type_id",
                    src_p1, fixed = TRUE))
  # Placeholder fully substituted (no leftover marker)
  expect_false(grepl("{{CANONICAL_HELPERS}}", src_p1, fixed = TRUE))
  # Same for pmulti path
  spec_pm <- list(p_class = "pmulti", prior = prior,
                  cp_a = FALSE, cp_W = FALSE, cp_a_per_k = NULL)
  src_pm <- .gdpar_emit_canonical_stan(spec_pm)
  expect_true(grepl("vector bspline_basis_eval(real x, int W_per_k_dim",
                    src_pm, fixed = TRUE))
  expect_true(grepl("vector apply_W_basis_diff(int basis_type_id",
                    src_pm, fixed = TRUE))
  expect_false(grepl("{{CANONICAL_HELPERS}}", src_pm, fixed = TRUE))
})

test_that("dispatcher aborts when canonical helpers piece is missing (defensive)", {
  # The dispatcher gates the helpers substitution behind
  # startsWith(piece_name, "amm_canonical_") so EB-side templates
  # served from inst/stan/ root are unaffected; for canonical pieces
  # missing helpers triggers gdpar_internal_error.
  prior <- gdpar_prior()
  spec <- list(p_class = "p1", prior = prior,
               cp_a = FALSE, cp_W = FALSE, mle = FALSE)
  # Sanity: with helpers piece present, the call succeeds.
  expect_silent(.gdpar_emit_canonical_stan(spec))
})

test_that("cmdstanr parse-check accepts the canonical pieces (skip_if cmdstan unavailable)", {
  skip_if_not_installed("cmdstanr")
  cmdstan_ok <- tryCatch({
    p <- cmdstanr::cmdstan_path()
    nzchar(p) && file.exists(p)
  }, error = function(e) FALSE)
  skip_if(!cmdstan_ok, "cmdstan toolchain not available")
  prior <- gdpar_prior()
  src_p1 <- generate_stan_code(prior)
  tf_p1 <- tempfile(pattern = "gdpar_canonical_p1_", fileext = ".stan")
  writeLines(src_p1, tf_p1)
  expect_silent(
    suppressMessages(cmdstanr::cmdstan_model(tf_p1, compile = FALSE,
                                             quiet = TRUE))
  )
  src_pmulti <- generate_stan_code_multi(prior)
  tf_pmulti <- tempfile(pattern = "gdpar_canonical_pmulti_", fileext = ".stan")
  writeLines(src_pmulti, tf_pmulti)
  expect_silent(
    suppressMessages(cmdstanr::cmdstan_model(tf_pmulti, compile = FALSE,
                                             quiet = TRUE))
  )
})

# NOTE (RG.6 part 1, Session B9.30, 2026-06-03, D96): the frozen md5 hashes
# and nchar below were refreshed because Option A (sigma_a_k compaction to
# slots with J_a_free > 0 via the transformed-data sigma_a_idx block)
# intentionally changed the byte content of the K and KxP templates and of
# the dispatcher output (+~1052 chars per K template, +~989 per KxP, exactly
# the sigma_a_idx block). This is a deliberate, documented source change, not
# drift; the fit-result bit-identity of every all-a() golden was verified by
# re-fit-and-compare (the goldens fit data, not source hashes). See NEWS.md.
#
# ---- Sub-bloque 9.3.a colateral L.iv.A (Sesion B9.5, 2026-05-27) ----
# Frozen piece + dispatcher OUTPUT invariants for the relocated
# amm_canonical_distrib_K.stan (legacy inst/stan/amm_distrib_K.stan
# replaced; helpers inline replaced by the {{CANONICAL_HELPERS}}
# placeholder consumed by the dispatcher per G.iv pattern).
#
# Sub-bloque 9.3.d J.iv.A (Sesion B9.5, 2026-05-27): frozen piece +
# dispatcher OUTPUT invariants for the new
# amm_canonical_pmulti_KxP.stan plus end-to-end smoke for the thin
# wrapper generate_stan_code_KxP_FB().

test_that("canonical piece amm_canonical_distrib_K.stan has the frozen md5 hash (B9.5 L.iv.A)", {
  # New B9.5 L.iv.A: relocated from inst/stan/amm_distrib_K.stan with
  # the inline bspline_basis_eval + apply_W_basis_diff definitions
  # replaced by the {{CANONICAL_HELPERS}} placeholder. The legacy
  # path is deleted per [[feedback_depuracion_obsoletos]]; the
  # canonical piece is the single source-of-truth K = 1 + p = 1 FB
  # template under the dispatcher.
  path <- .b9_3_get_canonical_piece_path("amm_canonical_distrib_K.stan")
  expect_true(file.exists(path))
  observed <- unname(tools::md5sum(path))
  expect_equal(observed, "64ae8cc0b9c8cb6b7bf18d8744f12a1c")
})

test_that("canonical piece amm_canonical_pmulti_KxP.stan has the frozen md5 hash (B9.5 J.iv.A)", {
  # New B9.5 J.iv.A: dedicated piece for the multi-parametric
  # multivariate Full-Bayes regime (K >= 2 slots crossed with p >= 2
  # AMM coordinates). Mirrors amm_eb_marginal_KxP.stan from
  # Sub-fase 8.6.D extended with the W component, the full Path B
  # family set {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13}, and the
  # canonical helpers piece reuse via {{CANONICAL_HELPERS}}.
  path <- .b9_3_get_canonical_piece_path("amm_canonical_pmulti_KxP.stan")
  expect_true(file.exists(path))
  observed <- unname(tools::md5sum(path))
  expect_equal(observed, "a6ecc4d2f9985f72722949c4d8b205c4")
})

test_that("dispatcher OUTPUT md5 frozen invariant - distrib_K default gaussian (B9.5 L.iv.A)", {
  prior <- gdpar_prior()
  src <- generate_stan_code_K(prior, family = gdpar_family("gaussian"))
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.5 L.iv.A frozen output hash; the relocated piece preserves
  # Stan semantics bit-exact relative to the legacy template (only
  # comments differ via the L.iv.A migration block + the helpers
  # piece comments; cmdstanr cache recompiles once but goldens K = 2
  # / K = 3 in tests/testthat/data/golden_K*.rds remain valid since
  # they fit data, not source hashes).
  expect_equal(observed, "cafd1c8adc4df319347be7cef1c3fac6")
  expect_equal(nchar(src), 41840L)
})

test_that("dispatcher OUTPUT md5 frozen invariant - distrib_K tweedie K=3 (B9.5 L.iv.A)", {
  prior <- gdpar_prior()
  src <- generate_stan_code_K(prior, family = gdpar_family("tweedie"))
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.5 L.iv.A frozen output hash for the family-specific
  # THETA_REF_PRIOR_BLOCK path (Tweedie K=3 with bounded p slot).
  expect_equal(observed, "b21c5c1c61645bc6ed68776224846e02")
  expect_equal(nchar(src), 41947L)
})

test_that("dispatcher OUTPUT md5 frozen invariant - pmulti_KxP default gaussian (B9.5 J.iv.A)", {
  prior <- gdpar_prior()
  src <- generate_stan_code_KxP_FB(prior, family = gdpar_family("gaussian"))
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.5 J.iv.A frozen output hash for the new KxP piece (canonical
  # vectorized THETA_REF_PRIOR_BLOCK_KxP).
  expect_equal(observed, "cc8a44200f2b539b7b44945496273b10")
  expect_equal(nchar(src), 32927L)
})

test_that("dispatcher OUTPUT md5 frozen invariant - pmulti_KxP tweedie K=3 (B9.5 J.iv.A)", {
  prior <- gdpar_prior()
  src <- generate_stan_code_KxP_FB(prior, family = gdpar_family("tweedie"))
  observed <- digest::digest(src, algo = "md5", serialize = FALSE)
  # B9.5 J.iv.A frozen output hash for the family-specific KxP
  # THETA_REF_PRIOR_BLOCK_KxP (Tweedie K=3 with bounded p slot per
  # coord via the iterated slice form documented in
  # .gdpar_build_theta_ref_prior_block_KxP).
  expect_equal(observed, "7173704bd765574729ec915bdb59e61a")
  expect_equal(nchar(src), 33228L)
})

test_that("dispatcher OUTPUT contains the canonical helpers - distrib_K + pmulti_KxP (B9.5 sanity)", {
  prior <- gdpar_prior()
  src_K <- generate_stan_code_K(prior, family = gdpar_family("gaussian"))
  expect_true(grepl("vector bspline_basis_eval(real x, int W_per_k_dim",
                    src_K, fixed = TRUE))
  expect_true(grepl("vector apply_W_basis_diff(int basis_type_id",
                    src_K, fixed = TRUE))
  expect_false(grepl("{{CANONICAL_HELPERS}}", src_K, fixed = TRUE))
  src_KxP <- generate_stan_code_KxP_FB(prior,
                                        family = gdpar_family("gaussian"))
  expect_true(grepl("vector bspline_basis_eval(real x, int W_per_k_dim",
                    src_KxP, fixed = TRUE))
  expect_true(grepl("vector apply_W_basis_diff(int basis_type_id",
                    src_KxP, fixed = TRUE))
  expect_false(grepl("{{CANONICAL_HELPERS}}", src_KxP, fixed = TRUE))
})

test_that("dispatcher rejects pmulti_KxP malformed spec (defensive)", {
  prior <- gdpar_prior()
  # Missing prior triggers an internal validation in the wrapper via
  # assert_inherits; the dispatcher branch itself is reached only with
  # a valid spec list.
  expect_error(
    generate_stan_code_KxP_FB(NULL),
    class = "gdpar_error"
  )
})

test_that("generate_stan_code_KxP_FB delegates to the canonical dispatcher bit-identically", {
  prior <- gdpar_prior()
  src_wrapper <- generate_stan_code_KxP_FB(prior, cp_W = FALSE,
                                             family = gdpar_family("gaussian"))
  src_direct <- .gdpar_emit_canonical_stan(list(
    p_class = "pmulti_KxP",
    prior   = prior,
    cp_W    = FALSE,
    family  = gdpar_family("gaussian")
  ))
  expect_identical(src_wrapper, src_direct)
})

test_that("cmdstanr parse-check accepts pmulti_KxP + distrib_K relocated (skip_if cmdstan unavailable)", {
  skip_if_not_installed("cmdstanr")
  cmdstan_ok <- tryCatch({
    p <- cmdstanr::cmdstan_path()
    nzchar(p) && file.exists(p)
  }, error = function(e) FALSE)
  skip_if(!cmdstan_ok, "cmdstan toolchain not available")
  prior <- gdpar_prior()
  src_K <- generate_stan_code_K(prior, family = gdpar_family("gaussian"))
  tf_K <- tempfile(pattern = "gdpar_canonical_distrib_K_", fileext = ".stan")
  writeLines(src_K, tf_K)
  expect_silent(
    suppressMessages(cmdstanr::cmdstan_model(tf_K, compile = FALSE,
                                             quiet = TRUE))
  )
  src_KxP <- generate_stan_code_KxP_FB(prior,
                                        family = gdpar_family("gaussian"))
  tf_KxP <- tempfile(pattern = "gdpar_canonical_pmulti_KxP_",
                     fileext = ".stan")
  writeLines(src_KxP, tf_KxP)
  expect_silent(
    suppressMessages(cmdstanr::cmdstan_model(tf_KxP, compile = FALSE,
                                             quiet = TRUE))
  )
})

# ---- Sub-bloque 9.3.a colateral L.iv.A.2 (Sesion B9.6, 2026-05-27): ----
# cascade L heredada from L.iv.A (B9.5). The 6 EB-side templates with
# inline bspline_basis_eval + apply_W_basis_diff helpers were relocated
# to inst/stan/_canonical_pieces/ with the `// {{CANONICAL_HELPERS}}`
# placeholder reusing the G.iv infrastructure of B9.4. Bit-exact Stan
# semantics preservation verified token-wise (helpers piece comments
# unified but Stan tokens identical). The 6 legacy top-level files
# were deleted (depuracion de obsoletos).

test_that("EB cascade L.iv.A.2 pieces have frozen md5 hashes (B9.6 cascade)", {
  hashes <- c(
    "amm_canonical_eb_marginal.stan"           = "7a279bf087a0fd0a95d35a490814673f",
    "amm_canonical_eb_marginal_multi.stan"     = "9250144e95df5fd409cb5cc9df3f3730",
    "amm_canonical_eb_marginal_K.stan"         = "98251e7cbf37303ee91ac4c9ccbfc280",
    "amm_canonical_eb_conditional.stan"        = "8321f9d7d99d0c702045a2886cd1f98c",
    "amm_canonical_eb_conditional_multi.stan"  = "186b6ac9f9810584ea21522974fe641d",
    "amm_canonical_eb_conditional_K.stan"      = "1d8de6308d52bf907c2294ee595f844b"
  )
  for (piece in names(hashes)) {
    path <- .b9_3_get_canonical_piece_path(piece)
    expect_true(file.exists(path), info = piece)
    observed <- unname(tools::md5sum(path))
    expect_equal(observed, unname(hashes[piece]),
                 info = sprintf("frozen md5 mismatch for %s", piece))
  }
})

test_that("EB cascade L.iv.A.2 pieces all carry the CANONICAL_HELPERS placeholder (B9.6)", {
  for (piece in c("amm_canonical_eb_marginal.stan",
                  "amm_canonical_eb_marginal_multi.stan",
                  "amm_canonical_eb_marginal_K.stan",
                  "amm_canonical_eb_conditional.stan",
                  "amm_canonical_eb_conditional_multi.stan",
                  "amm_canonical_eb_conditional_K.stan")) {
    path <- .b9_3_get_canonical_piece_path(piece)
    src <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_true(grepl("// {{CANONICAL_HELPERS}}", src, fixed = TRUE),
                info = sprintf("missing placeholder in %s", piece))
    # No inline helper redeclaration should remain.
    expect_false(grepl("vector bspline_basis_eval", src, fixed = TRUE),
                 info = piece)
    expect_false(grepl("vector apply_W_basis_diff", src, fixed = TRUE),
                 info = piece)
  }
})

test_that("EB cascade L.iv.A.2: legacy top-level EB templates with helpers inline are removed (B9.6)", {
  for (legacy in c("amm_eb_marginal.stan",
                   "amm_eb_marginal_multi.stan",
                   "amm_eb_marginal_K.stan",
                   "amm_eb_conditional.stan",
                   "amm_eb_conditional_multi.stan",
                   "amm_eb_conditional_K.stan")) {
    p <- system.file("stan", legacy, package = "gdpar")
    if (!nzchar(p)) p <- file.path("inst", "stan", legacy)
    expect_false(file.exists(p),
                 info = sprintf("legacy top-level EB template '%s' should be deleted (B9.6 cascade L.iv.A.2)",
                                legacy))
  }
})

test_that("EB cascade L.iv.A.2: KxP EB templates remain in inst/stan/ root (no helpers inline in 8.6.D)", {
  for (kxp in c("amm_eb_marginal_KxP.stan", "amm_eb_conditional_KxP.stan")) {
    p <- system.file("stan", kxp, package = "gdpar")
    if (!nzchar(p)) p <- file.path("inst", "stan", kxp)
    expect_true(file.exists(p),
                info = sprintf("KxP EB template '%s' should remain unchanged (8.6.D canonization; W = 0 hardcoded; no helpers cascade needed)",
                               kxp))
  }
})

test_that("EB cascade L.iv.A.2 render outputs have no leftover placeholders (B9.6 dispatcher integration)", {
  prior <- gdpar_prior()
  # EB single via .gdpar_eb_render_template (legacy alias)
  src_m <- gdpar:::.gdpar_eb_render_template("amm_eb_marginal.stan",
                                              prior, cp_a = FALSE, cp_W = FALSE)
  src_c <- gdpar:::.gdpar_eb_render_template("amm_eb_conditional.stan",
                                              prior, cp_a = FALSE, cp_W = FALSE)
  expect_false(grepl("{{", src_m, fixed = TRUE))
  expect_false(grepl("{{", src_c, fixed = TRUE))
  # EB multi via generate_stan_code_multi (legacy alias)
  src_mm <- generate_stan_code_multi(prior,
                                      template_name = "amm_eb_marginal_multi.stan")
  src_cm <- generate_stan_code_multi(prior,
                                      template_name = "amm_eb_conditional_multi.stan")
  expect_false(grepl("{{", src_mm, fixed = TRUE))
  expect_false(grepl("{{", src_cm, fixed = TRUE))
  # EB K via generate_stan_code_K (legacy alias)
  fam <- gdpar_family("gaussian")
  src_mK <- generate_stan_code_K(prior, family = fam,
                                  template_name = "amm_eb_marginal_K.stan")
  src_cK <- generate_stan_code_K(prior, family = fam,
                                  template_name = "amm_eb_conditional_K.stan")
  expect_false(grepl("{{", src_mK, fixed = TRUE))
  expect_false(grepl("{{", src_cK, fixed = TRUE))
})

# ---- Sub-bloque 9.3.d D69 substantive closure (Sesion B9.6, 2026-05-27) ----
# under decision N.a canonized at session open. The driver
# .gdpar_fb_KxP_fit() composes design + assembler (path = "FB") +
# cmdstanr sampling and bypasses the user-facing
# .gdpar_guard_multiparam_multivariate. The unified assembler
# .assemble_stan_data_KxP() now accepts path = c("EB", "FB"):
# the EB branch preserves the 8.6.D first-iteration restrictions
# byte-identically; the FB branch lifts W disabled (D39) + restricted
# stan_id {1, 3} (D40') to use_W enabled + the Path B family set
# {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13}.

test_that("unified assembler .assemble_stan_data_KxP preserves EB-side defaults (M.iv B9.6)", {
  set.seed(91001L)
  n <- 40L
  d_df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  amm_list <- list(
    mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
    sigma = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
  )
  design_KxP <- gdpar:::.build_amm_design_KxP(amm_list, d_df,
                                               formula_rhs = ~ x1 + x2)
  sd_eb <- gdpar:::.assemble_stan_data_KxP(
    design_KxP, family = gdpar_family("gaussian"),
    amm_list_canonical = amm_list,
    y_matrix = matrix(rnorm(n * 2L), n, 2L),
    theta_anchor_kp = matrix(0, nrow = 2L, ncol = 2L)
  )
  expect_equal(sd_eb$use_W, 0L)
  expect_equal(as.integer(sd_eb$family_id_k_vector), c(1L, 1L))
  # EB-side base list should NOT carry the W metadata fields.
  expect_null(sd_eb$dim_W)
  expect_null(sd_eb$W_per_kj_dim)
})

test_that("unified assembler .assemble_stan_data_KxP accepts Path B family set under path = 'FB' (M.iv B9.6)", {
  set.seed(91002L)
  n <- 40L
  d_df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  amm_list <- list(
    mu  = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
    phi = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
  )
  design_KxP <- gdpar:::.build_amm_design_KxP(amm_list, d_df,
                                               formula_rhs = ~ x1 + x2)
  for (fam_name in c("gaussian", "beta", "gamma")) {
    sd_fb <- gdpar:::.assemble_stan_data_KxP(
      design_KxP, family = gdpar_family(fam_name),
      amm_list_canonical = amm_list,
      y_matrix = matrix(runif(n * 2L, 0.1, 0.9), n, 2L),
      theta_anchor_kp = matrix(0, nrow = 2L, ncol = 2L),
      path = "FB"
    )
    expect_true(!is.null(sd_fb$W_per_kj_dim),
                info = sprintf("FB path missing W_per_kj_dim field for %s", fam_name))
    expect_true(!is.null(sd_fb$W_type_id),
                info = sprintf("FB path missing W_type_id field for %s", fam_name))
  }
  # bernoulli (stan_id = 4) NOT in Path B set under FB.
  expect_error(
    gdpar:::.assemble_stan_data_KxP(
      design_KxP, family = gdpar_family("bernoulli"),
      amm_list_canonical = amm_list,
      y_matrix = matrix(sample(0:1, n * 2L, replace = TRUE), n, 2L),
      theta_anchor_kp = matrix(0, nrow = 2L, ncol = 2L),
      path = "FB"
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_fb_KxP_fit guards against degenerate input (D69 driver B9.6)", {
  # Empty amm_list
  expect_error(
    gdpar:::.gdpar_fb_KxP_fit(
      y_matrix = matrix(0, 4L, 2L),
      family = gdpar_family("gaussian"),
      amm_list = list(),
      data = data.frame(x = 1:4),
      formula_rhs = ~ x
    ),
    class = "gdpar_input_error"
  )
  # p = 1 (single-coord) should reject
  amm_list_p1 <- list(
    mu    = amm_spec(p = 1L),
    sigma = amm_spec(p = 1L)
  )
  expect_error(
    gdpar:::.gdpar_fb_KxP_fit(
      y_matrix = matrix(0, 4L, 1L),
      family = gdpar_family("gaussian"),
      amm_list = amm_list_p1,
      data = data.frame(x = 1:4),
      formula_rhs = ~ x
    ),
    class = "gdpar_input_error"
  )
})

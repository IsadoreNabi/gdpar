# Tests for Sub-bloque 9.3.d (Bloque 9, Sesion B9.5 architectural
# canonization + Sesion B9.6 operational closure of deuda D69) under
# decisions I.iv + J.iv.A + M.iv + N.a. Mirror of
# tests/testthat/test-eb_goldens_8_6_D.R for the FB Path C bit-exact
# comparator template canonized in B9.5 via the new piece
# inst/stan/_canonical_pieces/amm_canonical_pmulti_KxP.stan plus the
# thin wrapper generate_stan_code_KxP_FB() that delegates to the
# canonical dispatcher .gdpar_emit_canonical_stan(), now coupled with
# the unified assembler .assemble_stan_data_KxP(path = "FB") and the
# fit driver .gdpar_fb_KxP_fit() canonized in B9.6.
#
# The actual cmdstanr sampling for the 4 reserved seeds 91001..91004
# is gated by env vars per the bootstrap script
# data-raw/bootstrap_fb_goldens_KxP.R; the RDS files appear in
# tests/testthat/data/ only when the bootstrap runs. Tests for the
# fitable RDS structure are now active (Section 4) once the
# bootstrap produces them. Sections 1-3 are codegen / bootstrap
# script / roster always-on or skip-on-absence.

.fb_KxP_roster_b9_5 <- list(
  list(name = "g_KxP_001", family = "gaussian",  K = 2L, p = 2L,
       seed = 91001L, basis_type = "polynomial"),
  list(name = "g_KxP_002", family = "beta",      K = 2L, p = 2L,
       seed = 91002L, basis_type = "polynomial"),
  list(name = "g_KxP_003", family = "gamma",     K = 2L, p = 3L,
       seed = 91003L, basis_type = "polynomial"),
  list(name = "g_KxP_004", family = "student_t", K = 3L, p = 2L,
       seed = 91004L, basis_type = "polynomial")
)

.fb_KxP_rds_path <- function(entry) {
  file.path("data",
            sprintf("golden_fb_KxP_%s_%s_K%d_p%d.rds",
                    entry$family, entry$basis_type, entry$K, entry$p))
}

.fb_KxP_first_golden_present <- function() {
  any(vapply(.fb_KxP_roster_b9_5,
             function(e) file.exists(.fb_KxP_rds_path(e)),
             logical(1L)))
}

# Section 1: codegen path tests (always run; do not depend on RDS) -----

test_that("FB KxP codegen emits the canonical piece for every roster entry (B9.5 J.iv.A)", {
  prior <- gdpar_prior()
  for (entry in .fb_KxP_roster_b9_5) {
    fam <- gdpar_family(entry$family)
    src <- generate_stan_code_KxP_FB(prior, family = fam)
    expect_true(nchar(src) > 5000L)
    expect_false(grepl("{{", src, fixed = TRUE))
    expect_true(grepl("array[J_groups, K] vector[p] theta_ref_kp",
                      src, fixed = TRUE))
    # Family-specific dispatch block present in the model{} section.
    fam_id <- as.integer(fam$stan_id)
    expect_true(
      grepl(sprintf("family_id_k_vector[1] == %d", fam_id),
            src, fixed = TRUE)
    )
  }
})

test_that("FB KxP codegen is deterministic per (family, seed) entry (B9.5)", {
  prior <- gdpar_prior()
  for (entry in .fb_KxP_roster_b9_5) {
    fam <- gdpar_family(entry$family)
    src_a <- generate_stan_code_KxP_FB(prior, family = fam)
    src_b <- generate_stan_code_KxP_FB(prior, family = fam)
    expect_identical(src_a, src_b)
  }
})

test_that("FB KxP reserved canonical seeds are inside the 91001..91099 range (B9.4 DESIGN)", {
  # DESIGN_9_3_D_PATH_C.md §4 reserved 91001..91099 for KxP goldens;
  # the K.c B9.5 roster occupies 91001..91004. The expanded roster
  # (K.c deferral to B9.6 / B9.x) MUST stay within 91005..91099.
  seeds <- vapply(.fb_KxP_roster_b9_5, function(e) e$seed, integer(1L))
  expect_true(all(seeds >= 91001L & seeds <= 91099L))
  expect_equal(length(unique(seeds)), length(seeds))
  expect_equal(seeds, c(91001L, 91002L, 91003L, 91004L))
})

# Section 2: bootstrap script presence + env-gating documentation ------

test_that("FB KxP bootstrap script exists and gates on canonical env vars (B9.5)", {
  # data-raw/ is excluded from the installed package by .Rbuildignore
  # default; locate via the source tree if available.
  script_path <- "../../data-raw/bootstrap_fb_goldens_KxP.R"
  if (!file.exists(script_path)) {
    script_path <- file.path("data-raw", "bootstrap_fb_goldens_KxP.R")
  }
  testthat::skip_if_not(file.exists(script_path),
                        "data-raw/bootstrap_fb_goldens_KxP.R not in source tree")
  src <- paste(readLines(script_path), collapse = "\n")
  expect_true(grepl("GDPAR_BOOTSTRAP_FB_GOLDENS", src, fixed = TRUE))
  expect_true(grepl("GDPAR_BOOTSTRAP_FB_GOLDENS_KxP", src, fixed = TRUE))
  expect_true(grepl("91001L", src, fixed = TRUE))
  expect_true(grepl("91002L", src, fixed = TRUE))
  expect_true(grepl("91003L", src, fixed = TRUE))
  expect_true(grepl("91004L", src, fixed = TRUE))
  expect_true(grepl("D69", src, fixed = TRUE))
})

# Section 3: golden RDS file presence (skip if not bootstrapped) -------

test_that("all FB KxP 9.3.d golden RDS files exist on disk (skip if not bootstrapped)", {
  testthat::skip_if_not(.fb_KxP_first_golden_present(),
                        "no FB KxP 9.3.d goldens present (bootstrap not yet run; deuda D69 pending B9.6)")
  for (entry in .fb_KxP_roster_b9_5) {
    expect_true(file.exists(.fb_KxP_rds_path(entry)),
                info = sprintf("missing golden RDS for %s", entry$name))
  }
})

test_that("FB KxP 9.3.d goldens carry the canonical metadata marker or fit object (skip if not bootstrapped)", {
  testthat::skip_if_not(.fb_KxP_first_golden_present(),
                        "no FB KxP 9.3.d goldens present")
  for (entry in .fb_KxP_roster_b9_5) {
    rds_path <- .fb_KxP_rds_path(entry)
    if (!file.exists(rds_path)) next
    payload <- readRDS(rds_path)
    if (is.list(payload) && !is.null(payload$status) &&
        payload$status == "metadata_only_pending_D69_B9_6") {
      # Legacy B9.5 stub regime (deprecated post-B9.6): metadata marker.
      expect_equal(payload$golden_id, entry$name)
      expect_equal(payload$seed, entry$seed)
      expect_equal(payload$family, entry$family)
      expect_equal(payload$K, entry$K)
      expect_equal(payload$p, entry$p)
    } else {
      # B9.6 regime: a fitable snapshot list with draws_matrix +
      # canonical metadata (status = "fitable").
      expect_true(is.list(payload))
      expect_equal(payload$golden_id, entry$name)
      expect_equal(payload$seed, entry$seed)
      expect_equal(payload$family, entry$family)
      expect_equal(payload$K, entry$K)
      expect_equal(payload$p, entry$p)
      expect_equal(payload$sub_phase, "9.3.d")
      expect_equal(payload$status, "fitable")
      expect_true(!is.null(payload$draws_matrix))
      expect_true(!is.null(payload$stan_data))
    }
  }
})

# Section 3b: Option A sigma_a_k compaction guard (RG.6, D96) ------------

test_that("Option A (D96): g_KxP_004 (nu a = ~ 1) compacts sigma_a_k below K", {
  testthat::skip_if_not(.fb_KxP_first_golden_present(),
                        "no FB KxP 9.3.d goldens present")
  st <- Filter(function(e) e$family == "student_t", .fb_KxP_roster_b9_5)
  testthat::skip_if(length(st) == 0L, "no student_t KxP roster entry")
  entry <- st[[1L]]
  rds_path <- .fb_KxP_rds_path(entry)
  testthat::skip_if_not(file.exists(rds_path), "g_KxP_004 golden not present")
  payload <- readRDS(rds_path)
  if (is.list(payload) && !is.null(payload$status) &&
      payload$status == "metadata_only_pending_D69_B9_6") {
    testthat::skip("g_KxP_004 still a metadata stub")
  }
  cn <- colnames(payload$draws_matrix)
  sak <- grep("^sigma_a_k\\[", cn, value = TRUE)
  # The nu slot declares an intercept-only a() (a = ~ 1) -> it carries no
  # free a coefficient -> its sigma_a_k scale is a flat direction
  # (non-identification) that Option A (J_a_free > 0 criterion, D96)
  # removes. With K = 3 slots and only mu / sigma carrying free a
  # coefficients, the compacted sigma_a_k has length 2; the dropped
  # third scale is the faithful marginalization validated at re-bootstrap.
  expect_length(sak, 2L)
  expect_setequal(sak, c("sigma_a_k[1]", "sigma_a_k[2]"))
  expect_false("sigma_a_k[3]" %in% cn)
})

# Section 4: Tier 1 fitable structure (B9.6 D69 closure) -----------------

test_that("FB KxP 9.3.d fitable goldens have non-empty draws_matrix with theta_ref_kp (B9.6 N.a)", {
  testthat::skip_if_not(.fb_KxP_first_golden_present(),
                        "no FB KxP 9.3.d goldens present")
  for (entry in .fb_KxP_roster_b9_5) {
    rds_path <- .fb_KxP_rds_path(entry)
    if (!file.exists(rds_path)) next
    payload <- readRDS(rds_path)
    if (is.list(payload) && !is.null(payload$status) &&
        payload$status == "metadata_only_pending_D69_B9_6") next
    dm <- payload$draws_matrix
    expect_true(NROW(dm) > 0L,
                info = sprintf("empty draws_matrix in %s", entry$name))
    # The KxP piece declares theta_ref_kp as
    # array[J_groups, K] vector[p]; cmdstanr serializes it with the
    # canonical "theta_ref_kp[g,k,j]" name pattern. Verify at least
    # one element exists.
    cn <- colnames(dm)
    has_theta <- any(grepl("^theta_ref_kp\\[", cn))
    expect_true(has_theta,
                info = sprintf("no theta_ref_kp column in %s draws", entry$name))
  }
})

test_that("FB KxP 9.3.d fitable goldens stan_data carries the FB extension fields (B9.6 M.iv)", {
  testthat::skip_if_not(.fb_KxP_first_golden_present(),
                        "no FB KxP 9.3.d goldens present")
  for (entry in .fb_KxP_roster_b9_5) {
    rds_path <- .fb_KxP_rds_path(entry)
    if (!file.exists(rds_path)) next
    payload <- readRDS(rds_path)
    if (is.list(payload) && !is.null(payload$status) &&
        payload$status == "metadata_only_pending_D69_B9_6") next
    sd <- payload$stan_data
    # FB path output must carry the M.iv extension fields.
    for (f in c("use_W", "dim_W", "d", "W_per_kj_dim", "W_type_id",
                "W_n_knots_full", "W_knots_full", "W_degree")) {
      expect_true(f %in% names(sd),
                  info = sprintf("missing FB extension field %s in %s",
                                 f, entry$name))
    }
    # Sanity: K + p match roster.
    expect_equal(as.integer(sd$K), as.integer(entry$K),
                 info = entry$name)
    expect_equal(as.integer(sd$p), as.integer(entry$p),
                 info = entry$name)
  }
})

test_that("manifest CSV contains sub_phase = '9.3.d' rows when bootstrapped", {
  manifest_path <- gdpar:::.gdpar_golden_manifest_path()
  testthat::skip_if_not(file.exists(manifest_path),
                        "golden_manifest.csv not present")
  testthat::skip_if_not(.fb_KxP_first_golden_present(),
                        "no FB KxP 9.3.d goldens present")
  manifest <- gdpar:::.gdpar_golden_manifest_read(manifest_path)
  rows <- manifest[manifest$sub_phase == "9.3.d", , drop = FALSE]
  expect_true(nrow(rows) >= 1L)
  # Seeds present in manifest must be a subset of the reserved
  # 91001..91099 range.
  seeds_in_manifest <- suppressWarnings(as.integer(rows$seed))
  seeds_in_manifest <- seeds_in_manifest[!is.na(seeds_in_manifest)]
  if (length(seeds_in_manifest) > 0L) {
    expect_true(all(seeds_in_manifest >= 91001L &
                    seeds_in_manifest <= 91099L))
  }
})

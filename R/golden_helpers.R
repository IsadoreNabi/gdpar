# Sub-phase 8.3.9 (2026-05-22): unified bootstrap helpers + manifest
# CSV scaffolding for K=2/K=3 goldens (decisions D1=(1D),
# D2=(2B), E3=(E3.A) of 8.3.9 scoping).
#
# The manifest CSV lives at tests/testthat/data/golden_manifest.csv
# and tracks reproducibility metadata for every golden snapshot:
# family, K, basis_type, p, sub_phase, bootstrap_date,
# cmdstan_version, R_version, DHARMa_version, fit_code_hash, seed,
# n_chains, n_iter_warmup, n_iter_sampling. The fit_code_hash is the
# SHA256 of the deparse of the fit-code closure that produced the
# snapshot, so the auditor can detect drift in the generator.
#
# Actual bootstrap of the 14 configurations is gated by env vars per
# golden (one BOOTSTRAP flag per family x K x basis_type) so the
# whole suite stays fast in CI; an instance with NOT_CRAN=true and
# GDPAR_GOLDEN_BOOTSTRAP_ALL=1 can rebuild every snapshot in one go.

#' Internal: path to the golden manifest CSV.
#'
#' @keywords internal
#' @noRd
.gdpar_golden_manifest_path <- function(root_dir = NULL) {
  if (is.null(root_dir)) {
    candidates <- c("tests/testthat/data", "data", ".")
    for (c in candidates) {
      if (dir.exists(c)) {
        return(file.path(c, "golden_manifest.csv"))
      }
    }
    return("golden_manifest.csv")
  }
  file.path(root_dir, "golden_manifest.csv")
}

#' Internal: canonical column order of the manifest.
#'
#' @keywords internal
#' @noRd
.gdpar_golden_manifest_columns <- function() {
  c("family", "K", "basis_type", "p", "sub_phase", "bootstrap_date",
    "cmdstan_version", "R_version", "DHARMa_version", "fit_code_hash",
    "seed", "n_chains", "n_iter_warmup", "n_iter_sampling",
    "rds_path", "n_obs")
}

#' Internal: initialize the manifest CSV with the canonical header if
#' the file does not exist; otherwise no-op.
#'
#' @keywords internal
#' @noRd
.gdpar_golden_manifest_init <- function(path = NULL) {
  if (is.null(path)) path <- .gdpar_golden_manifest_path()
  if (file.exists(path)) return(invisible(path))
  header <- .gdpar_golden_manifest_columns()
  empty <- as.data.frame(
    matrix(character(0), nrow = 0L, ncol = length(header),
           dimnames = list(NULL, header)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(empty, path, row.names = FALSE)
  invisible(path)
}

#' Internal: read the manifest CSV (returns a data.frame; empty if
#' the file does not exist).
#'
#' @keywords internal
#' @noRd
.gdpar_golden_manifest_read <- function(path = NULL) {
  if (is.null(path)) path <- .gdpar_golden_manifest_path()
  if (!file.exists(path)) {
    cols <- .gdpar_golden_manifest_columns()
    return(as.data.frame(
      matrix(character(0), nrow = 0L, ncol = length(cols),
             dimnames = list(NULL, cols)),
      stringsAsFactors = FALSE
    ))
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Internal: append one row to the manifest CSV (idempotent on the
#' tuple family x K x basis_type x p: replaces any existing row with
#' the same key).
#'
#' @keywords internal
#' @noRd
.gdpar_golden_manifest_add <- function(row_data, path = NULL) {
  if (is.null(path)) path <- .gdpar_golden_manifest_path()
  .gdpar_golden_manifest_init(path)
  current <- .gdpar_golden_manifest_read(path)
  required <- c("family", "K", "basis_type", "p")
  missing_keys <- setdiff(required, names(row_data))
  if (length(missing_keys) > 0L) {
    gdpar_abort(
      sprintf(
        "Golden manifest row is missing required keys: %s.",
        paste(sQuote(missing_keys), collapse = ", ")
      ),
      class = "gdpar_internal_error"
    )
  }
  cols <- .gdpar_golden_manifest_columns()
  row_df <- as.data.frame(
    matrix(NA_character_, nrow = 1L, ncol = length(cols),
           dimnames = list(NULL, cols)),
    stringsAsFactors = FALSE
  )
  for (k in names(row_data)) {
    if (k %in% cols) row_df[[k]][1L] <- as.character(row_data[[k]])
  }
  if (nrow(current) > 0L) {
    key_old <- paste(current$family, current$K, current$basis_type,
                     current$p, sep = "|")
    key_new <- paste(row_df$family[1L], row_df$K[1L],
                     row_df$basis_type[1L], row_df$p[1L], sep = "|")
    keep <- key_old != key_new
    current <- current[keep, , drop = FALSE]
  }
  combined <- rbind(current, row_df)
  utils::write.csv(combined, path, row.names = FALSE)
  invisible(path)
}

#' Internal: compute SHA256 hash of the fit-code closure (deparse +
#' digest) as the reproducibility identifier for a golden.
#'
#' Falls back to a string representation when the \pkg{digest}
#' package is unavailable.
#'
#' @keywords internal
#' @noRd
.gdpar_golden_fit_code_hash <- function(fit_code_fn) {
  txt <- paste(deparse(fit_code_fn), collapse = "\n")
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(txt, algo = "sha256"))
  }
  paste0("nohash_", as.character(nchar(txt)))
}

#' Internal: cmdstan version string (or NA if cmdstanr not installed).
#'
#' @keywords internal
#' @noRd
.gdpar_cmdstan_version <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    return(NA_character_)
  }
  v <- tryCatch(cmdstanr::cmdstan_version(),
                error = function(e) NA_character_)
  as.character(v)
}

#' Internal: DHARMa version string (or NA).
#'
#' @keywords internal
#' @noRd
.gdpar_dharma_version <- function() {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion("DHARMa"))
}

#' Internal: canonical filename for a golden RDS keyed by
#' family x K x basis_type x p (matches the convention of the
#' three frozen goldens of 8.3.4 / 8.3.5a / 8.3.3 Unit 3).
#'
#' @keywords internal
#' @noRd
.gdpar_golden_rds_name <- function(family, K, basis_type = "polynomial",
                                   p = 1L) {
  base <- sprintf("golden_K%d_%s_K%d", K, family, K)
  if (basis_type != "polynomial") {
    base <- paste(base, basis_type, sep = "_")
  }
  if (p > 1L) {
    base <- sprintf("%s_p%d", base, p)
  }
  paste0(base, ".rds")
}

#' Internal: canonical filename for an EB golden RDS keyed by
#' the (family, configuration, status) triple of Sub-phase 8.6.B.
#'
#' The convention is parallel to \code{.gdpar_golden_rds_name}: the
#' filename embeds the family name plus the regime label (a string
#' identifying the configuration, e.g. \code{"poly"}, \code{"bspline"},
#' \code{"guard"}). Examples:
#'   \itemize{
#'     \item \code{golden_eb_gaussian_poly.rds}
#'     \item \code{golden_eb_gaussian_bspline.rds}
#'     \item \code{golden_eb_zip_guard.rds}
#'   }
#'
#' @keywords internal
#' @noRd
.gdpar_golden_eb_rds_name <- function(family, regime = "poly") {
  paste0("golden_eb_", family, "_", regime, ".rds")
}

#' Internal: canonical 17-configuration roster for Sub-phase 8.6.B
#' (Charter Section 3.2 / parcial handoff Section 5.1 + D33 of the closure session).
#'
#' 5 fitable (status = "fitable") + 12 guard-rejected
#' (status = "guard"). Each entry carries: \code{family} (canonical
#' name), \code{regime} (filename suffix), \code{status} (fitable /
#' guard), \code{guard_class} (NA for fitable; the canonical
#' \code{gdpar_unsupported_feature_error} for the guarded ones),
#' \code{rationale} (one-line description of the guard rationale or
#' the fit-time canonical regime).
#'
#' @keywords internal
#' @noRd
.gdpar_golden_roster_8_6_B <- function() {
  list(
    # ---- Fitable (5) -----------------------------------------------------
    list(family = "gaussian",       regime = "poly",   status = "fitable",
         guard_class = NA_character_,
         rationale = "K=1, polynomial W, identity link"),
    list(family = "poisson",        regime = "poly",   status = "fitable",
         guard_class = NA_character_,
         rationale = "K=1, polynomial W, log link"),
    list(family = "neg_binomial_2", regime = "poly",   status = "fitable",
         guard_class = NA_character_,
         rationale = "K=1 with phi as global dispersion, log link"),
    list(family = "bernoulli",      regime = "poly",   status = "fitable",
         guard_class = NA_character_,
         rationale = "K=1, polynomial W, logit link"),
    list(family = "gaussian",       regime = "bspline", status = "fitable",
         guard_class = NA_character_,
         rationale = "K=1, B-spline W (Sub-phase 8.3.8 path)"),
    # ---- Guard-rejected, stan_id D33 (5) ---------------------------------
    list(family = "beta",          regime = "guard",  status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=5 outside EB-supported set (D33)"),
    list(family = "gamma",         regime = "guard",  status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=6 outside EB-supported set (D33)"),
    list(family = "student_t",     regime = "guard",  status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=8 outside EB-supported set (D33)"),
    list(family = "tweedie",       regime = "guard",  status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=9 outside EB-supported set (D33)"),
    list(family = "lognormal_loc_scale", regime = "guard",
         status = "guard",
         guard_class = "gdpar_input_error",
         rationale = "not in gdpar_family() registry (K>=2 only)"),
    # ---- Guard-rejected, min_K (4) ---------------------------------------
    list(family = "zip",                  regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "min_K>1 (mixture; deferred to 8.6.C)"),
    list(family = "zinb",                 regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "min_K>1 (mixture; deferred to 8.6.C)"),
    list(family = "hurdle_poisson",       regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "min_K>1 (hurdle; deferred to 8.6.C)"),
    list(family = "hurdle_neg_binomial_2", regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "min_K>1 (hurdle; deferred to 8.6.C)"),
    # ---- Guard-rejected, heterogeneous amm K>1 (3) -----------------------
    list(family = "het_beta_gamma_K2",    regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "heterogeneous K=2 amm list; deferred to 8.6.C"),
    list(family = "het_tweedie_K3",       regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "heterogeneous K=3 amm list; deferred to 8.6.C"),
    list(family = "het_student_t_K3",     regime = "guard",
         status = "guard",
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "heterogeneous K=3 amm list; deferred to 8.6.C")
  )
}

#' Internal: canonical 17-configuration roster for Sub-phase 8.6.C
#' (Charter Section 3.3, apertura Section 4.4, parcial / Session 11
#' D34 canonization). Covers the two multivariate paths habilitados
#' por D34:
#' \itemize{
#'   \item Path A (K = 1, p > 1): 6 fitable (Gaussian / Poisson /
#'     Negative Binomial / Bernoulli on amm_eb_marginal_multi.stan;
#'     stan_id in {1, 2, 3, 4}) + 4 guard (Beta / Gamma / Student-t /
#'     Tweedie on the Path A check_stan_id_for_path because Path A
#'     supported set inherits amm_distrib_multi.stan {1, 2, 3, 4}).
#'   \item Path B (K > 1, p = 1): 7 fitable spanning Beta K=2, Gamma
#'     K=2, Student-t K=3, Tweedie K=3, Lognormal-loc-scale K=2,
#'     ZIP K=2 and heterogeneous Gauss+Beta K=2.
#' }
#'
#' Each entry carries: \code{family} (canonical name; for hetero K
#' uses the joint label \code{"het_gauss_beta_K2"}), \code{regime}
#' (filename suffix, \code{"polyP2"} / \code{"polyP3"} /
#' \code{"bsplineP2"} / \code{"K2"} / \code{"K3"} /
#' \code{"pathA_p2_guard"}), \code{status} (\code{"fitable"} /
#' \code{"guard"}), \code{path} (\code{"A"} or \code{"B"}),
#' \code{K} / \code{p} (integers), \code{guard_class} (NA for
#' fitable; canonical error class for guarded entries),
#' \code{rationale}.
#'
#' Filename pattern: \code{golden_eb_<family>_path<A|B>_<regime>.rds}.
#'
#' @keywords internal
#' @noRd
.gdpar_golden_roster_8_6_C <- function() {
  list(
    # ---- Path A K=1 p>1 fitable (6) --------------------------------------
    list(family = "gaussian", regime = "polyP2",
         status = "fitable", path = "A", K = 1L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path A K=1 p=2 polynomial Gaussian (stan_id 1)"),
    list(family = "gaussian", regime = "bsplineP2",
         status = "fitable", path = "A", K = 1L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path A K=1 p=2 B-spline Gaussian (stan_id 1)"),
    list(family = "gaussian", regime = "polyP3",
         status = "fitable", path = "A", K = 1L, p = 3L,
         guard_class = NA_character_,
         rationale = "Path A K=1 p=3 polynomial Gaussian (stan_id 1)"),
    list(family = "poisson", regime = "polyP2",
         status = "fitable", path = "A", K = 1L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path A K=1 p=2 polynomial Poisson (stan_id 2)"),
    list(family = "neg_binomial_2", regime = "polyP2",
         status = "fitable", path = "A", K = 1L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path A K=1 p=2 polynomial NB2 (stan_id 3)"),
    list(family = "bernoulli", regime = "polyP2",
         status = "fitable", path = "A", K = 1L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path A K=1 p=2 polynomial Bernoulli (stan_id 4)"),
    # ---- Path A K=1 p>1 guard (4) ----------------------------------------
    list(family = "beta", regime = "polyP2_guard",
         status = "guard", path = "A", K = 1L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=5 outside Path A supported set (K=1)"),
    list(family = "gamma", regime = "polyP2_guard",
         status = "guard", path = "A", K = 1L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=6 outside Path A supported set (K=1)"),
    list(family = "student_t", regime = "polyP2_guard",
         status = "guard", path = "A", K = 1L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=8 outside Path A supported set (K=1)"),
    list(family = "tweedie", regime = "polyP2_guard",
         status = "guard", path = "A", K = 1L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=9 outside Path A supported set (K=1)"),
    # ---- Path B K>1 p=1 fitable (7) --------------------------------------
    list(family = "beta", regime = "K2",
         status = "fitable", path = "B", K = 2L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=2 polynomial Beta (stan_id 5)"),
    list(family = "gamma", regime = "K2",
         status = "fitable", path = "B", K = 2L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=2 polynomial Gamma (stan_id 6)"),
    list(family = "student_t", regime = "K3",
         status = "fitable", path = "B", K = 3L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=3 polynomial Student-t (stan_id 8)"),
    list(family = "tweedie", regime = "K3",
         status = "fitable", path = "B", K = 3L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=3 polynomial Tweedie (stan_id 9)"),
    list(family = "lognormal_loc_scale", regime = "K2",
         status = "fitable", path = "B", K = 2L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=2 lognormal-loc-scale (stan_id 7)"),
    list(family = "zip", regime = "K2",
         status = "fitable", path = "B", K = 2L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=2 ZIP mixture (stan_id 10)"),
    list(family = "het_gauss_beta_K2", regime = "K2",
         status = "fitable", path = "B", K = 2L, p = 1L,
         guard_class = NA_character_,
         rationale = "Path B K=2 heterogeneous Gauss + Beta (named-list family)")
  )
}

#' Internal: canonical filename for an EB 8.6.C golden RDS keyed by
#' the (family, path, regime) triple. Parallel to
#' \code{.gdpar_golden_eb_rds_name} but embedding the Path A / Path B
#' tag in the basename so the 8.6.B and 8.6.C namespaces do not
#' collide on disk.
#'
#' Examples:
#'   \itemize{
#'     \item \code{golden_eb_gaussian_pathA_polyP2.rds}
#'     \item \code{golden_eb_beta_pathB_K2.rds}
#'     \item \code{golden_eb_het_gauss_beta_K2_pathB_K2.rds}
#'   }
#'
#' @keywords internal
#' @noRd
.gdpar_golden_eb_rds_name_8_6_C <- function(family, path, regime) {
  sprintf("golden_eb_%s_path%s_%s.rds", family, path, regime)
}

#' Internal: canonical roster for Sub-phase 8.6.D (Path C K>1 + p>1)
#'
#' Sub-phase 8.6.D (decisions D36 = alpha + D37 = (i) Session 13a +
#' D38'' = (h) + D39 + D40' + D41 Session 13b + D43 = (a) Session 13c).
#' Initial iteration coverage: Path C \code{stan_id in {1, 3}} (Gaussian
#' K=2 + NB K=2) crossed with p \code{in {2, 3}}. Other Path B families
#' \code{{5, 6, 7, 8, 9, 10, 11, 12, 13}} are guard-rejected entries
#' documenting the deferred status with the explicit numerical caveat
#' of opening Section 6.1.
#'
#' Total: 4 fitable (Gaussian / NB x p in {2, 3}) + 4 metadata-only
#' guard (Beta / Gamma / Student-t / ZIP at p = 2 with the deferred
#' message). Heterogeneous-p across slots is registered as a follow-on
#' debt against Block 9.x and not represented in the initial roster.
#'
#' Each entry carries: \code{family} (canonical name), \code{regime}
#' (filename suffix, e.g. \code{"polyP2"} / \code{"polyP3"} /
#' \code{"polyP2_guard"}), \code{status} (\code{"fitable"} /
#' \code{"guard"}), \code{path} (\code{"C"}), \code{K} / \code{p}
#' (integers), \code{guard_class} (NA for fitable; canonical error
#' class for guarded entries), \code{rationale}.
#'
#' Filename pattern: \code{golden_eb_<family>_pathC_<regime>.rds} via
#' \code{.gdpar_golden_eb_rds_name_8_6_D}.
#'
#' @keywords internal
#' @noRd
.gdpar_golden_roster_8_6_D <- function() {
  list(
    # ---- Path C fitable (4) ----------------------------------------------
    list(family = "gaussian", regime = "polyP2",
         status = "fitable", path = "C", K = 2L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path C K=2 p=2 polynomial Gaussian K=2 (stan_id 1)"),
    list(family = "gaussian", regime = "polyP3",
         status = "fitable", path = "C", K = 2L, p = 3L,
         guard_class = NA_character_,
         rationale = "Path C K=2 p=3 polynomial Gaussian K=2 (stan_id 1)"),
    list(family = "neg_binomial_2", regime = "polyP2",
         status = "fitable", path = "C", K = 2L, p = 2L,
         guard_class = NA_character_,
         rationale = "Path C K=2 p=2 polynomial NB K=2 (stan_id 3)"),
    list(family = "neg_binomial_2", regime = "polyP3",
         status = "fitable", path = "C", K = 2L, p = 3L,
         guard_class = NA_character_,
         rationale = "Path C K=2 p=3 polynomial NB K=2 (stan_id 3)"),
    # ---- Path C guard (4) -----------------------------------------------
    list(family = "beta", regime = "polyP2_guard",
         status = "guard", path = "C", K = 2L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=5 deferred under D40' (numerical caveat 6.1)"),
    list(family = "gamma", regime = "polyP2_guard",
         status = "guard", path = "C", K = 2L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=6 deferred under D40' (numerical caveat 6.1)"),
    list(family = "student_t", regime = "polyP2_guard",
         status = "guard", path = "C", K = 3L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=8 deferred under D40' (numerical caveat 6.1)"),
    list(family = "zip", regime = "polyP2_guard",
         status = "guard", path = "C", K = 2L, p = 2L,
         guard_class = "gdpar_unsupported_feature_error",
         rationale = "stan_id=10 deferred under D40' (numerical caveat 6.1)")
  )
}

#' Internal: canonical filename for an EB 8.6.D golden RDS keyed by
#' the (family, path, regime) triple. Mirrors
#' \code{.gdpar_golden_eb_rds_name_8_6_C} with the path tag hardcoded
#' to "C" so the 8.6.B / 8.6.C / 8.6.D namespaces do not collide on
#' disk.
#'
#' Examples:
#'   \itemize{
#'     \item \code{golden_eb_gaussian_pathC_polyP2.rds}
#'     \item \code{golden_eb_neg_binomial_2_pathC_polyP3.rds}
#'     \item \code{golden_eb_beta_pathC_polyP2_guard.rds}
#'   }
#'
#' @keywords internal
#' @noRd
.gdpar_golden_eb_rds_name_8_6_D <- function(family, regime) {
  sprintf("golden_eb_%s_pathC_%s.rds", family, regime)
}

#' Internal: canonical 14-configuration roster for 8.3.9 (Decision
#' D1 = (1D), 13 atrasados + K=1+p=2 polynomial baseline).
#'
#' @keywords internal
#' @noRd
.gdpar_golden_roster_8_3_9 <- function() {
  list(
    # 8.3.4 atrasados: custom_lognormal K=2 (gaussian + beta + gamma
    # already exist in disk)
    list(family = "lognormal_loc_scale", K = 2L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.4"),
    # 8.3.5a Student-t K=3
    list(family = "student_t", K = 3L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.5a"),
    # 8.3.5b Tweedie K=3
    list(family = "tweedie", K = 3L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.5b"),
    # 8.3.6 mixtures (4)
    list(family = "zip", K = 2L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.6"),
    list(family = "zinb", K = 3L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.6"),
    list(family = "hurdle_poisson", K = 2L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.6"),
    list(family = "hurdle_neg_binomial_2", K = 3L,
         basis_type = "polynomial", p = 1L, sub_phase = "8.3.6"),
    # 8.3.7 heterogeneous K=2
    list(family = "het_gauss_beta", K = 2L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.7"),
    list(family = "het_gauss_gamma", K = 2L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.7"),
    list(family = "het_nb_beta", K = 2L, basis_type = "polynomial",
         p = 1L, sub_phase = "8.3.7"),
    # 8.3.8 bspline (3)
    list(family = "gaussian", K = 1L, basis_type = "bspline",
         p = 1L, sub_phase = "8.3.8"),
    list(family = "het_gauss_beta", K = 2L, basis_type = "bspline",
         p = 1L, sub_phase = "8.3.8"),
    list(family = "gaussian", K = 1L, basis_type = "bspline",
         p = 2L, sub_phase = "8.3.8"),
    # 8.3.9 baseline: K=1 + p=2 polynomial (filling the last gap)
    list(family = "gaussian", K = 1L, basis_type = "polynomial",
         p = 2L, sub_phase = "8.3.9")
  )
}

#' Internal: count the number of individual-scope parameters of a family
#'
#' Returns the integer K = number of \code{param_specs} entries whose
#' \code{scope} is per-observation or per-group. Population auxiliaries
#' (e.g. \code{sigma}, \code{phi} in the Block 8.0 default K = 1) do not
#' count. Block 8 design decision 2a-iii (hybrid scope) and decision 1C
#' (family as list of param_specs).
#'
#' Used by the guard against K > 1 + p > 1 input combinations: in Block
#' 8.0 the multi-parametric extension is exposed only for p = 1; the
#' multivariate path requires K = 1 until Block 8.1 lifts the guard.
#'
#' For backward compatibility with family objects that pre-date the
#' \code{param_specs} field (Block 8.0 sub-phase a refactor), the
#' function falls back to K = 1 when \code{param_specs} is absent.
#'
#' @param family Either a \code{gdpar_family} or \code{gdpar_family_multi}.
#' @return Integer scalar K >= 1.
#' @keywords internal
#' @noRd
.gdpar_n_params_individual <- function(family) {
  if (inherits(family, "gdpar_family_multi")) {
    base <- family$families[[1L]]
  } else {
    base <- family
  }
  if (is.null(base$param_specs)) {
    return(1L)
  }
  individual_scopes <- c("per_observation", "per_group")
  K <- sum(vapply(base$param_specs,
                  function(s) s$scope %in% individual_scopes,
                  logical(1L)))
  as.integer(max(K, 1L))
}

#' Internal: guard against K > 1 with p > 1 in Block 8.0
#'
#' Block 8.0 of the multi-parametric extension exposes K > 1 only for
#' p = 1; K > 1 with p > 1 is the combinatoric regime queued for Block
#' 8.1 under decision (M-iii'). The guard aborts with an informative
#' error when a user-supplied family + AMM specification crosses into
#' that regime.
#'
#' @param family A \code{gdpar_family} or \code{gdpar_family_multi}.
#' @param p Integer scalar with the AMM dimension.
#' @return Invisibly NULL when the guard passes; aborts otherwise.
#' @keywords internal
#' @noRd
.gdpar_guard_multiparam_multivariate <- function(family, p) {
  K <- .gdpar_n_params_individual(family)
  if (K > 1L && p > 1L) {
    gdpar_abort(
      paste(
        "K > 1 (multi-parametric family) combined with p > 1 (multivariate",
        "AMM) is deferred to Block 8.1. In Block 8.0 the multi-parametric",
        "extension is exposed only for p = 1; the multivariate AMM remains",
        "single-parameter per coordinate (K = 1). Use either K > 1 with",
        "p = 1, or K = 1 with p > 1."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(K = K, p = p)
    )
  }
  invisible(NULL)
}

#' Internal: guard against routing a family below its minimum required K
#'
#' Some built-in families are wired exclusively under the K-individual
#' distributional regression path of \code{amm_distrib_K.stan} because
#' their structural auxiliary parameters (dispersion, degrees of
#' freedom, ...) are part of the likelihood definition. The legacy
#' K = 1 templates (\code{amm_main.stan},
#' \code{amm_distrib_multi.stan}) do not branch on their
#' \code{stan_id}; routing such a family below its minimum slot count
#' would yield a Stan model with no likelihood for those branches. The
#' guard aborts before that with an informative pointer to the
#' elevation via \code{gdpar_bf()} or the equivalent named-list of
#' \code{amm_spec} objects.
#'
#' Minimum K per stan_id (post sub-phase 8.3.5b, 2026-05-21):
#' \itemize{
#'   \item stan_id 5 (Beta): minimum K = 2 (mu + phi).
#'   \item stan_id 6 (Gamma): minimum K = 2 (mu + shape).
#'   \item stan_id 7 (lognormal_loc_scale, custom-K registry): minimum K = 2 (mu + sigma).
#'   \item stan_id 8 (Student-t): minimum K = 3 (mu + sigma + nu).
#'   \item stan_id 9 (Tweedie): minimum K = 3 (mu + phi + p).
#'   \item stan_id 10 (ZIP): minimum K = 2 (mu + pi).
#'   \item stan_id 11 (ZINB): minimum K = 3 (mu + phi + pi).
#'   \item stan_id 12 (hurdle_poisson): minimum K = 2 (mu + pi).
#'   \item stan_id 13 (hurdle_neg_binomial_2): minimum K = 3
#'     (mu + phi + pi).
#' }
#' Families with K = 1 wiring in the legacy templates (stan_id 1, 2,
#' 3, 4) bypass the guard.
#'
#' @param family A \code{gdpar_family} or \code{gdpar_family_multi}.
#' @return Invisibly NULL when the guard passes; aborts otherwise.
#' @keywords internal
#' @noRd
.gdpar_guard_K_below_family_min <- function(family) {
  if (inherits(family, "gdpar_family_multi")) {
    base <- family$families[[1L]]
  } else {
    base <- family
  }
  if (is.null(base$stan_id) || is.na(base$stan_id)) {
    return(invisible(NULL))
  }
  min_K <- switch(as.character(base$stan_id),
    "5" = 2L,
    "6" = 2L,
    "7" = 2L,
    "8" = 3L,
    "9" = 3L,
    "10" = 2L,
    "11" = 3L,
    "12" = 2L,
    "13" = 3L,
    NULL
  )
  if (is.null(min_K)) {
    return(invisible(NULL))
  }
  K <- .gdpar_n_params_individual(family)
  if (K >= min_K) {
    return(invisible(NULL))
  }
  slot_names <- if (!is.null(base$param_specs)) {
    vapply(base$param_specs, function(s) s$name, character(1L))
  } else {
    c("location", "dispersion")
  }
  aux_slots <- if (length(slot_names) >= 2L) slot_names[-1L] else "dispersion"
  aux_summary <- paste(sQuote(aux_slots), collapse = ", ")
  example_first_aux <- aux_slots[1L]
  gdpar_abort(
    sprintf(
      paste(
        "Family '%s' requires K = %d K-individual slots under the",
        "distributional regression path of amm_distrib_K.stan (received",
        "K = %d). Elevate the auxiliary slot(s) %s via the formula",
        "interface, e.g.\n",
        "    gdpar(gdpar_bf(y ~ a(x1), %s ~ a(x2), ...),\n",
        "          family = gdpar_family('%s'), ...)\n",
        "or the equivalent named-list of amm_spec objects. Routings",
        "with K < %d are queued for a future sub-phase."
      ),
      base$name, min_K, K, aux_summary, example_first_aux,
      base$name, min_K
    ),
    class = "gdpar_unsupported_feature_error",
    data = list(family_name = base$name, stan_id = base$stan_id,
                min_K = min_K, received_K = K,
                aux_slots = aux_slots)
  )
}

#' Internal: canonical dispatcher for the K=1 FB Stan templates
#'
#' Single source-of-truth dispatcher for the (p, K=1) FB Stan templates
#' introduced in Sub-sub-fase 9.3.a (Bloque 9, Sesion B9.3, 2026-05-27)
#' under the canonized decision B.iv lateral (codegen R-side as
#' canonical source of truth). Reads the canonical piece selected by
#' \code{spec$p_class} from \code{inst/stan/_canonical_pieces/} and
#' applies the appropriate placeholder substitution. Replaces the
#' duplicated template-read plus substitute logic that previously lived
#' inline in \code{generate_stan_code()} (p = 1) and
#' \code{generate_stan_code_multi()} (p >= 1). The pieces in
#' \code{_canonical_pieces/} are frozen byte-identical copies of the
#' legacy \code{amm_main.stan} and \code{amm_distrib_multi.stan}; the
#' bit-exactness of the substituted Stan source is preserved by
#' construction (no golden re-bootstrap required) since the dispatcher
#' emits the same substituted strings as the legacy paths for the same
#' inputs.
#'
#' Templates whose name starts with \code{"amm_canonical_"} are served
#' from \code{inst/stan/_canonical_pieces/}; other names (EB-side
#' templates pending cascade in B9.5+) are served from
#' \code{inst/stan/} root.
#'
#' Sub-sub-fase 9.3.a colateral (Sesion B9.4, 2026-05-27) under
#' decision G.iv lateral: when the body piece contains the
#' \code{// \{\{CANONICAL_HELPERS\}\}} placeholder, the dispatcher
#' reads \code{amm_canonical_helpers.stan} (a Stan source FRAGMENT
#' shipped in \code{_canonical_pieces/}) and substitutes it into the
#' placeholder via \code{gsub()} before applying the other
#' parametrization placeholders. The helpers piece holds the canonical
#' definitions of \code{bspline_basis_eval} and \code{apply_W_basis_diff}
#' (deuda f of 9.3.a, closed); the body pieces only retain the
#' placeholder inside their \code{functions \{ \}} block.
#'
#' @param spec Named list with the canonical codegen specification.
#'   Required fields by \code{p_class}:
#'   \describe{
#'     \item{\code{p_class = "p1"}}{Single-coordinate template
#'       (\code{amm_canonical_p1.stan}). Required: \code{prior},
#'       \code{cp_a}, \code{cp_W}, \code{mle}.}
#'     \item{\code{p_class = "pmulti"}}{Multi-coordinate template
#'       (\code{amm_canonical_pmulti.stan} or an EB-side template that
#'       shares the placeholder set bit-for-bit). Required:
#'       \code{prior}, \code{cp_a}, \code{cp_W}, \code{cp_a_per_k}
#'       (may be \code{NULL} to broadcast \code{cp_a}),
#'       \code{template_name} (defaults to
#'       \code{"amm_canonical_pmulti.stan"}).}
#'     \item{\code{p_class = "pmulti_KxP"}}{Multi-parametric
#'       multivariate template (\code{amm_canonical_pmulti_KxP.stan})
#'       canonized in Sub-bloque 9.3.d (Bloque 9, Sesion B9.5,
#'       2026-05-27) under decision I.iv lateral + J.iv.A piece
#'       arquitectonica (DESIGN_9_3_D_PATH_C.md sub-decision 3.2).
#'       Covers K >= 2 K-individual slots crossed with p >= 2 AMM
#'       coordinates; reuses the canonical helpers piece via
#'       \code{// \{\{CANONICAL_HELPERS\}\}}. Required: \code{prior},
#'       \code{cp_a} (NCP uniform default; per-slot per-coord CP /
#'       mixed deferred to B9.6+), \code{cp_W}, \code{family} (NULL
#'       defaults to the canonical vectorized THETA_REF_PRIOR_BLOCK;
#'       Tweedie K=3 uses the bounded p-slot slice form).}
#'   }
#' @return Character scalar containing the substituted Stan source
#'   ready for \code{cmdstanr::cmdstan_model}.
#' @keywords internal
#' @noRd
.gdpar_emit_canonical_stan <- function(spec) {
  if (!is.list(spec) || is.null(spec$p_class) ||
      !spec$p_class %in% c("p1", "pmulti", "distrib_K",
                            "pmulti_KxP")) {
    gdpar_abort(
      sprintf(
        paste(
          "Invalid codegen spec: p_class must be 'p1', 'pmulti',",
          "'distrib_K' or 'pmulti_KxP', got '%s'."
        ),
        format(spec$p_class)
      ),
      class = "gdpar_internal_error"
    )
  }
  if (spec$p_class == "p1") {
    piece_name <- "amm_canonical_p1.stan"
  } else if (spec$p_class == "pmulti_KxP") {
    piece_name <- "amm_canonical_pmulti_KxP.stan"
  } else if (spec$p_class == "distrib_K") {
    piece_name <- if (is.null(spec$template_name)) {
      "amm_canonical_distrib_K.stan"
    } else {
      spec$template_name
    }
  } else {
    piece_name <- if (is.null(spec$template_name)) {
      "amm_canonical_pmulti.stan"
    } else {
      spec$template_name
    }
  }
  if (startsWith(piece_name, "amm_canonical_")) {
    template_path <- system.file("stan", "_canonical_pieces", piece_name,
                                 package = "gdpar")
    if (!nzchar(template_path)) {
      template_path <- file.path("inst", "stan", "_canonical_pieces",
                                 piece_name)
    }
  } else {
    template_path <- system.file("stan", piece_name, package = "gdpar")
    if (!nzchar(template_path)) {
      template_path <- file.path("inst", "stan", piece_name)
    }
  }
  if (!file.exists(template_path)) {
    gdpar_abort(
      sprintf("Stan piece '%s' not found in the installed package.",
              piece_name),
      class = "gdpar_internal_error"
    )
  }
  src <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
  # Sub-sub-fase 9.3.a colateral (Bloque 9, Sesion B9.4, 2026-05-27)
  # under decision G.iv lateral: inject the canonical Stan helpers
  # (bspline_basis_eval + apply_W_basis_diff) from the dedicated piece
  # `amm_canonical_helpers.stan` via R-side textual substitution at the
  # `// {{CANONICAL_HELPERS}}` placeholder inside the body piece's
  # `functions { }` block. Closes deuda f of 9.3.a (helpers duplication
  # across pieces) without cmdstan `#include` semantics: the helpers
  # piece is a Stan source FRAGMENT (no surrounding `functions { }`)
  # and is concatenated R-side into the body piece pre-substitution.
  # The Stan semantics of the dispatcher output are preserved bit-exact
  # by construction; comments in the helpers piece are unified across
  # body pieces (the cmdstanr cache hash recompiles once on first call
  # post-B9.4 but the model is semantically identical).
  if (grepl("// {{CANONICAL_HELPERS}}", src, fixed = TRUE)) {
    helpers_path <- system.file("stan", "_canonical_pieces",
                                "amm_canonical_helpers.stan",
                                package = "gdpar")
    if (!nzchar(helpers_path)) {
      helpers_path <- file.path("inst", "stan", "_canonical_pieces",
                                "amm_canonical_helpers.stan")
    }
    if (!file.exists(helpers_path)) {
      gdpar_abort(
        "Canonical helpers piece 'amm_canonical_helpers.stan' not found.",
        class = "gdpar_internal_error"
      )
    }
    helpers_src <- paste(readLines(helpers_path, warn = FALSE),
                         collapse = "\n")
    src <- gsub("// {{CANONICAL_HELPERS}}", helpers_src, src, fixed = TRUE)
  }
  prior <- spec$prior
  cp_W <- isTRUE(spec$cp_W)
  if (spec$p_class == "p1") {
    cp_a <- isTRUE(spec$cp_a)
    if (isTRUE(spec$mle)) {
      cp_a <- TRUE
      cp_W <- TRUE
    }
    replacements <- list(
      "{{PRIOR_THETA_REF}}"       = prior$theta_ref,
      "{{PRIOR_SIGMA_THETA_REF}}" = prior$sigma_theta_ref,
      "{{PRIOR_SIGMA_A}}"         = prior$sigma_a,
      "{{PRIOR_SIGMA_B}}"         = prior$sigma_b,
      "{{PRIOR_SIGMA_W}}"         = prior$sigma_W,
      "{{PRIOR_SIGMA_Y}}"         = prior$sigma_y,
      "{{PRIOR_PHI}}"             = prior$phi,
      "{{A_SCALE}}"               = if (cp_a) "" else " * sigma_a[1]",
      "{{A_PRIOR}}"               = if (cp_a) "normal(0, sigma_a[1])"
                                           else "normal(0, 1)",
      "{{W_SCALE}}"               = if (cp_W) "" else " * sigma_W[1]",
      "{{W_PRIOR}}"               = if (cp_W) "normal(0, sigma_W[1])"
                                           else "normal(0, 1)"
    )
  } else if (spec$p_class == "distrib_K") {
    # Sub-bloque 9.3.a colateral L.iv.A (Bloque 9, Sesion B9.5,
    # 2026-05-27): the distrib_K branch consumes the same placeholder
    # set as the legacy generate_stan_code_K, including
    # THETA_REF_PRIOR_BLOCK (per-family anchor block via
    # .gdpar_build_theta_ref_prior_block) and DATA_CP_A_PER_K_DECL /
    # TP_A_BLOCK / MODEL_A_BLOCK (per-slot a-component via
    # generate_a_blocks_K with the cp_a_per_K logical vector). The
    # template_name override supports the EB-side K templates
    # (amm_eb_marginal_K.stan, amm_eb_conditional_K.stan) per
    # D34 of Sub-fase 8.6.C; those EB templates still live in
    # inst/stan/ root pending the cascade in B9.6 (deuda L heredada
    # bajo L.iv.A symmetry).
    cp_a_per_K <- spec$cp_a_per_K
    if (is.null(cp_a_per_K)) {
      cp_a_per_K <- as.logical(spec$cp_a)[1L]
    } else if (!is.logical(cp_a_per_K) || length(cp_a_per_K) < 1L ||
               any(is.na(cp_a_per_K))) {
      gdpar_abort(
        "Argument 'cp_a_per_K' must be a non-empty logical vector with no NA.",
        class = "gdpar_input_error",
        data = list(received = cp_a_per_K)
      )
    }
    blocks <- generate_a_blocks_K(cp_a_per_K)
    family <- spec$family
    theta_ref_prior_block <-
      .gdpar_build_theta_ref_prior_block(family)
    replacements <- list(
      "{{THETA_REF_PRIOR_BLOCK}}"  = theta_ref_prior_block,
      "{{PRIOR_THETA_REF}}"        = prior$theta_ref,
      "{{PRIOR_SIGMA_THETA_REF}}"  = prior$sigma_theta_ref,
      "{{PRIOR_SIGMA_A}}"          = prior$sigma_a,
      "{{PRIOR_SIGMA_B}}"          = prior$sigma_b,
      "{{PRIOR_SIGMA_W}}"          = prior$sigma_W,
      "{{PRIOR_SIGMA_Y}}"          = prior$sigma_y,
      "{{PRIOR_PHI}}"              = prior$phi,
      "{{DATA_CP_A_PER_K_DECL}}"   = blocks$data_decl,
      "{{TP_A_BLOCK}}"             = blocks$tp_block,
      "{{MODEL_A_BLOCK}}"          = blocks$model_block,
      "{{W_SCALE}}"                = if (cp_W) "" else " * sigma_W[1]",
      "{{W_PRIOR}}"                = if (cp_W) "normal(0, sigma_W[1])"
                                            else "normal(0, 1)"
    )
  } else if (spec$p_class == "pmulti_KxP") {
    # Sub-bloque 9.3.d (Bloque 9, Sesion B9.5, 2026-05-27) under
    # decision I.iv + J.iv.A: the KxP piece uses inline NCP-only a/b
    # blocks (no DATA_CP_A_PER_K_DECL / TP_A_BLOCK / MODEL_A_BLOCK
    # substitutions for B9.5 atomic scope; per-slot per-coord CP /
    # mixed variants deferred to B9.6+ as deuda emergente). The
    # placeholder set is the union of the canonical priors and the
    # THETA_REF_PRIOR_BLOCK shared with the K piece.
    family <- spec$family
    theta_ref_prior_block <-
      .gdpar_build_theta_ref_prior_block_KxP(family)
    replacements <- list(
      "{{THETA_REF_PRIOR_BLOCK}}"  = theta_ref_prior_block,
      "{{PRIOR_THETA_REF}}"        = prior$theta_ref,
      "{{PRIOR_SIGMA_THETA_REF}}"  = prior$sigma_theta_ref,
      "{{PRIOR_SIGMA_A}}"          = prior$sigma_a,
      "{{PRIOR_SIGMA_B}}"          = prior$sigma_b,
      "{{PRIOR_SIGMA_W}}"          = prior$sigma_W,
      "{{PRIOR_SIGMA_Y}}"          = prior$sigma_y,
      "{{PRIOR_PHI}}"              = prior$phi,
      "{{W_SCALE}}"                = if (cp_W) "" else " * sigma_W[1]",
      "{{W_PRIOR}}"                = if (cp_W) "normal(0, sigma_W[1])"
                                            else "normal(0, 1)"
    )
  } else {
    cp_a_per_k <- spec$cp_a_per_k
    if (is.null(cp_a_per_k)) {
      cp_a_per_k <- as.logical(spec$cp_a)[1L]
    } else if (!is.logical(cp_a_per_k) || length(cp_a_per_k) < 1L ||
               any(is.na(cp_a_per_k))) {
      gdpar_abort(
        "Argument 'cp_a_per_k' must be a non-empty logical vector with no NA.",
        class = "gdpar_input_error",
        data = list(received = cp_a_per_k)
      )
    }
    blocks <- generate_a_blocks_multi(cp_a_per_k)
    replacements <- list(
      "{{PRIOR_THETA_REF}}"        = prior$theta_ref,
      "{{PRIOR_SIGMA_THETA_REF}}"  = prior$sigma_theta_ref,
      "{{PRIOR_SIGMA_A}}"          = prior$sigma_a,
      "{{PRIOR_SIGMA_B}}"          = prior$sigma_b,
      "{{PRIOR_SIGMA_W}}"          = prior$sigma_W,
      "{{PRIOR_SIGMA_Y}}"          = prior$sigma_y,
      "{{PRIOR_PHI}}"              = prior$phi,
      "{{DATA_CP_A_PER_K_DECL}}"   = blocks$data_decl,
      "{{TP_A_BLOCK}}"             = blocks$tp_block,
      "{{MODEL_A_BLOCK}}"          = blocks$model_block,
      "{{W_SCALE}}"                = if (cp_W) "" else " * sigma_W[1]",
      "{{W_PRIOR}}"                = if (cp_W) "normal(0, sigma_W[1])"
                                            else "normal(0, 1)"
    )
  }
  for (k in names(replacements)) {
    src <- gsub(k, replacements[[k]], src, fixed = TRUE)
  }
  if (grepl("{{", src, fixed = TRUE)) {
    leftover <- regmatches(
      src, regexpr("\\{\\{[A-Za-z0-9_]+\\}\\}", src)
    )
    gdpar_abort(
      sprintf(
        "Unsubstituted placeholder remains in generated canonical Stan code: %s",
        leftover
      ),
      class = "gdpar_internal_error"
    )
  }
  if (spec$p_class == "p1" && isTRUE(spec$mle)) {
    src <- .strip_priors_block(src)
  }
  src
}

#' Internal: read the Path 1 Stan template and substitute prior placeholders
#'
#' Loads the static Stan template shipped with the package and replaces
#' the prior placeholders (e.g., \code{\{\{PRIOR_THETA_REF\}\}}) with
#' the user-supplied prior strings, plus the parametrization
#' placeholders that select between centered and non-centered forms
#' for the additive and modulating components. The substituted result
#' is the Stan model fitted by \code{cmdstanr}.
#'
#' Caching is delegated to \code{cmdstanr}, which hashes the source
#' string and reuses the compiled binary across calls when the source
#' is unchanged. CP and NCP variants therefore compile separately and
#' cache independently.
#'
#' @param prior An object of class \code{gdpar_prior}.
#' @param cp_a Logical scalar. If \code{TRUE}, the additive component
#'   \code{a} is sampled in the centered parametrization
#'   (\code{a_raw[j] ~ normal(0, sigma_a[1])} and
#'   \code{a_coef[j] = a_raw[j]}). If \code{FALSE} (default), the
#'   non-centered parametrization is used
#'   (\code{a_raw[j] ~ normal(0, 1)} and
#'   \code{a_coef[j] = a_raw[j] * sigma_a[1]}). The flag is determined
#'   automatically by the pre-flight diagnostic in \code{gdpar()}.
#' @param cp_W Logical scalar. If \code{TRUE}, the modulating
#'   component \code{W} is sampled in the centered parametrization
#'   (\code{W_raw[k, j] ~ normal(0, sigma_W[1])} with no scaling at
#'   use site). If \code{FALSE} (default), the non-centered
#'   parametrization is used (\code{W_raw[k, j] ~ normal(0, 1)} and
#'   the scaling \code{W_raw[k, j] * sigma_W[1]} is applied inside
#'   the polynomial-basis loop in transformed parameters). Determined
#'   automatically by the pre-flight diagnostic in \code{gdpar()}.
#' @param mle Logical scalar. If \code{TRUE}, emit the maximum-likelihood
#'   variant: every line between the markers \code{// BEGIN PRIORS} and
#'   \code{// END PRIORS} in the template is stripped, and \code{cp_a},
#'   \code{cp_W} are forced to \code{TRUE} (centered parametrization is
#'   the only one that yields an identified likelihood for \code{a_coef},
#'   \code{W_coef} when the hyper-priors are absent; without priors,
#'   the NCP product \code{a_raw[j] * sigma_a[1]} would be unidentified).
#'   Intended for \code{gdpar_bvm_check()} only; not used in the
#'   standard inference path. Defaults to \code{FALSE}.
#' @return A character scalar containing the Stan model source ready
#'   for \code{cmdstanr::cmdstan_model(write_stan_file(...))}.
#' @keywords internal
#' @noRd
generate_stan_code <- function(prior, cp_a = FALSE, cp_W = FALSE,
                                mle = FALSE) {
  assert_inherits(prior, "gdpar_prior", "prior")
  .gdpar_emit_canonical_stan(list(
    p_class = "p1",
    prior = prior,
    cp_a = cp_a,
    cp_W = cp_W,
    mle = mle
  ))
}

#' Internal: strip the prior block between BEGIN PRIORS / END PRIORS markers
#'
#' Removes every line from \code{// BEGIN PRIORS} through \code{// END PRIORS}
#' inclusive. Used by \code{generate_stan_code(mle = TRUE)} to emit the
#' MLE variant of the model consumed by \code{gdpar_bvm_check}. Aborts
#' if the markers are absent or malformed (template invariant violated).
#'
#' @param src Character scalar with the substituted Stan source.
#' @return Character scalar without the prior block.
#' @keywords internal
#' @noRd
.strip_priors_block <- function(src) {
  lines <- strsplit(src, "\n", fixed = TRUE)[[1L]]
  begin_idx <- grep("//\\s*BEGIN PRIORS", lines)
  end_idx   <- grep("//\\s*END PRIORS",   lines)
  if (length(begin_idx) != 1L || length(end_idx) != 1L ||
      begin_idx >= end_idx) {
    gdpar_abort(
      paste(
        "Stan template invariant violated:",
        "expected exactly one '// BEGIN PRIORS' and one '// END PRIORS'",
        "marker with BEGIN before END."
      ),
      class = "gdpar_internal_error"
    )
  }
  kept <- lines[-(begin_idx:end_idx)]
  paste(kept, collapse = "\n")
}

#' Internal: write generated Stan code to a temporary file
#'
#' CRAN policy requires that any file the package writes lives in
#' \code{tempdir()}. This helper writes the Stan model source to a
#' temporary file and returns the path, suitable for
#' \code{cmdstanr::cmdstan_model}. The file name embeds a hash of the
#' source so identical sources reuse the same path within a session,
#' enabling \code{cmdstanr}'s caching.
#'
#' @param src Character scalar with the Stan source code.
#' @return Character scalar with the path to the written file.
#' @keywords internal
#' @noRd
write_stan_to_tempfile <- function(src) {
  hash <- substr(
    paste0(format(as.numeric(charToRaw(src))), collapse = ""),
    1L, 16L
  )
  hash <- gsub("[^A-Za-z0-9]", "0", hash)
  path <- file.path(tempdir(), paste0("gdpar_amm_", hash, ".stan"))
  writeLines(src, path)
  path
}

#' Internal: assemble the data block passed to Stan
#'
#' Builds the named list consumed by \code{cmdstanr::sample}'s
#' \code{data} argument from the AMM design matrices, the family
#' object, the outcome vector and the anchor value.
#'
#' @param design List returned by \code{build_amm_design}.
#' @param family An object of class \code{gdpar_family}.
#' @param amm An object of class \code{amm_spec}.
#' @param y Numeric or integer vector of outcomes (length n).
#' @param theta_anchor Numeric scalar with the anchor value.
#' @return Named list with the entries expected by the Stan template.
#' @keywords internal
#' @noRd
assemble_stan_data <- function(design, family, amm, y, theta_anchor,
                                group_id = NULL) {
  p_value <- if (!is.null(amm$p)) as.integer(amm$p) else 1L
  .gdpar_guard_multiparam_multivariate(family, p_value)
  .gdpar_guard_K_below_family_min(family)
  if (inherits(family, "gdpar_family_multi") ||
      (!is.null(amm$p) && amm$p > 1L)) {
    return(.assemble_stan_data_multi(design, family, amm, y, theta_anchor,
                                      group_id = group_id))
  }
  n <- length(y)
  grp <- .resolve_group_id(group_id, n)
  use_a <- as.integer(!is.null(amm$a))
  use_b <- as.integer(!is.null(amm$b))
  use_W <- as.integer(!is.null(amm$W))
  J_a <- ncol(design$Z_a)
  J_b <- ncol(design$Z_b)
  d <- ncol(design$X)
  dim_W <- if (use_W == 1) amm$W$dim else 0L
  W_stan <- .gdpar_resolve_W_stan_data(
    W = amm$W, use_W = use_W, y_vec = y, family = family
  )
  needs_real <- family$name == "gaussian"
  needs_int <- family$name %in% c("poisson", "neg_binomial_2",
                                   "bernoulli")
  y_real <- if (needs_real) as.double(y) else rep(0.0, n)
  y_int <- if (needs_int) as.integer(y) else rep(0L, n)
  Z_a_mat <- if (J_a > 0) design$Z_a else matrix(0, n, 0)
  Z_b_mat <- if (J_b > 0) design$Z_b else matrix(0, n, 0)
  X_mat <- if (d > 0) design$X else matrix(0, n, 0)
  use_dispersion_y <- as.integer(family$name == "gaussian")
  use_dispersion_phi <- as.integer(family$name == "neg_binomial_2")
  list(
    n = n,
    family_id = family$stan_id,
    use_a = use_a,
    use_b = use_b,
    use_W = use_W,
    J_a = J_a,
    J_b = J_b,
    dim_W = as.integer(dim_W),
    d = as.integer(d),
    Z_a = Z_a_mat,
    Z_b = Z_b_mat,
    X = X_mat,
    y_real = y_real,
    y_int = y_int,
    theta_anchor = as.double(theta_anchor),
    use_dispersion_y = use_dispersion_y,
    use_dispersion_phi = use_dispersion_phi,
    W_type_id      = W_stan$W_type_id,
    W_n_knots_full = W_stan$W_n_knots_full,
    W_knots_full   = W_stan$W_knots_full,
    W_degree       = W_stan$W_degree,
    use_groups = grp$use_groups,
    J_groups   = grp$J_groups,
    group_id   = grp$group_id
  )
}

#' Internal: resolve the group_id argument into Stan-side data fields
#'
#' Returns a list with use_groups (0 or 1), J_groups (>= 1) and an
#' integer vector group_id of length n with values in 1..J_groups.
#' When group_id_arg is NULL the result is the backward-compatible
#' single-group regime; otherwise the input is coerced to a factor
#' to assign contiguous integer levels.
#'
#' @keywords internal
#' @noRd
.resolve_group_id <- function(group_id_arg, n) {
  if (is.null(group_id_arg)) {
    return(list(use_groups = 0L, J_groups = 1L,
                group_id = rep(1L, n)))
  }
  if (length(group_id_arg) != n) {
    gdpar_abort(
      sprintf(
        paste(
          "Argument 'group_id' must have length n = %d; received %d."
        ),
        n, length(group_id_arg)
      ),
      class = "gdpar_input_error",
      data = list(n = n, received_length = length(group_id_arg))
    )
  }
  if (any(is.na(group_id_arg))) {
    gdpar_abort(
      "Argument 'group_id' must not contain NA values.",
      class = "gdpar_input_error"
    )
  }
  fac <- as.factor(group_id_arg)
  list(use_groups = 1L,
       J_groups   = nlevels(fac),
       group_id   = as.integer(fac))
}

#' Internal: assemble the data block passed to the multivariate Stan template
#'
#' Multivariate path of \code{assemble_stan_data}. Consumes the ragged
#' design returned by \code{.build_amm_design_multi} and the
#' multivariate family object, and produces the data list expected by
#' the \code{amm_distrib_multi.stan} template: per-coordinate column
#' counts (\code{J_a_per_k}, \code{J_b_per_k}), padded design matrices
#' (\code{Z_a}, \code{Z_b} as \code{array[p] matrix[n, J_max]}), the
#' shared modulating block (\code{X}, \code{dim_W}, \code{W_per_k_dim}),
#' and per-coordinate outcomes (\code{y_real}, \code{y_int} as
#' \code{matrix[n, p]} / \code{array[n, p] int}).
#'
#' Stan does not support natively ragged arrays. The padding strategy
#' fills missing columns of each \code{Z_*[k]} with zeros up to
#' \code{J_*_max}; the Stan inner loops iterate over
#' \code{1:J_*_per_k[k]} so the padded entries never enter the
#' likelihood. The flat-pack of the per-coordinate coefficient vectors
#' (\code{a_raw}, \code{c_b_raw}) happens inside the Stan template via
#' offsets computed in \code{transformed data}; this assembler only
#' needs to expose the per-coordinate column counts.
#'
#' @param design List returned by \code{.build_amm_design_multi}.
#' @param family An object of class \code{gdpar_family_multi}.
#' @param amm An object of class \code{amm_spec} with \code{p > 1L}.
#' @param y Numeric matrix of dimensions \code{n} by \code{p} carrying
#'   the outcomes. Each column \eqn{k} feeds the marginal family
#'   \eqn{D_k(y_{ik} \mid \theta_i[k])}.
#' @param theta_anchor Either a scalar (broadcast to length \code{p})
#'   or a numeric vector of length \code{p} with the anchor values
#'   per coordinate.
#' @return Named list with the entries expected by
#'   \code{amm_distrib_multi.stan}.
#' @keywords internal
#' @noRd
.assemble_stan_data_multi <- function(design, family, amm, y, theta_anchor,
                                       group_id = NULL) {
  assert_inherits(amm, "amm_spec", "amm")
  if (is.null(amm$p) || amm$p < 1L) {
    gdpar_abort(
      "Internal error: amm$p must be set and positive for the multivariate assembler.",
      class = "gdpar_internal_error"
    )
  }
  p <- amm$p
  if (!inherits(family, "gdpar_family_multi")) {
    gdpar_abort(
      paste(
        "When amm$p > 1L, 'family' must be of class gdpar_family_multi.",
        "Use gdpar_family_multi() to construct it."
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(family))
    )
  }
  if (family$p != p) {
    gdpar_abort(
      sprintf(
        "Family p = %d does not match amm$p = %d. They must coincide.",
        family$p, p
      ),
      class = "gdpar_input_error",
      data = list(family_p = family$p, amm_p = p)
    )
  }
  if (!is.matrix(y) && !is.array(y)) {
    gdpar_abort(
      sprintf(
        paste(
          "When amm$p > 1L (multivariate path), 'y' must be a matrix of",
          "dimensions n by p = %d."
        ),
        p
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(y))
    )
  }
  y_mat <- if (is.matrix(y)) y else as.matrix(y)
  if (ncol(y_mat) != p) {
    gdpar_abort(
      sprintf(
        paste(
          "'y' must have %d columns (one per coordinate of theta_i);",
          "received %d."
        ),
        p, ncol(y_mat)
      ),
      class = "gdpar_input_error",
      data = list(expected_cols = p, received_cols = ncol(y_mat))
    )
  }
  n <- nrow(y_mat)
  any_a <- any(vapply(amm$dims, function(d) !is.null(d$a), logical(1L)))
  any_b <- any(vapply(amm$dims, function(d) !is.null(d$b), logical(1L)))
  use_a <- as.integer(any_a)
  use_b <- as.integer(any_b)
  use_W <- as.integer(!is.null(amm$W))
  J_a_per_k <- vapply(design$Z_a_list, ncol, integer(1L))
  J_b_per_k <- vapply(design$Z_b_list, ncol, integer(1L))
  J_a_max <- if (length(J_a_per_k) > 0L) max(J_a_per_k) else 0L
  J_b_max <- if (length(J_b_per_k) > 0L) max(J_b_per_k) else 0L
  pad_to <- function(z, target_cols, n_rows) {
    if (target_cols == 0L) {
      return(matrix(0, nrow = n_rows, ncol = 0L))
    }
    if (ncol(z) >= target_cols) {
      return(z)
    }
    cbind(z, matrix(0, nrow = nrow(z), ncol = target_cols - ncol(z)))
  }
  Z_a_padded <- lapply(design$Z_a_list, pad_to,
                       target_cols = J_a_max, n_rows = n)
  Z_b_padded <- lapply(design$Z_b_list, pad_to,
                       target_cols = J_b_max, n_rows = n)
  d <- ncol(design$X)
  if (use_W == 1L) {
    dim_W <- amm$W$dim
    if (is.null(dim_W) || !is.numeric(dim_W) || length(dim_W) != 1L ||
        dim_W <= 0L) {
      gdpar_abort(
        "Internal error: amm$W$dim must be a positive scalar for the multivariate assembler.",
        class = "gdpar_internal_error",
        data = list(dim_W = dim_W)
      )
    }
    if (dim_W %% p != 0L) {
      gdpar_abort(
        sprintf(
          paste(
            "dim_W = %d must be a positive multiple of p = %d for",
            "separable bases."
          ),
          as.integer(dim_W), p
        ),
        class = "gdpar_internal_error",
        data = list(dim_W = dim_W, p = p)
      )
    }
    W_per_k_dim <- as.integer(dim_W / p)
  } else {
    dim_W <- 0L
    W_per_k_dim <- 0L
  }
  W_stan <- .gdpar_resolve_W_stan_data(
    W = if (use_W == 1L) amm$W else NULL,
    use_W = use_W,
    y_vec = as.numeric(y_mat),
    family = family
  )
  needs_real <- family$name == "gaussian"
  needs_int <- family$name %in% c("poisson", "neg_binomial_2", "bernoulli")
  y_real <- if (needs_real) {
    storage.mode(y_mat) <- "double"
    y_mat
  } else {
    matrix(0, nrow = n, ncol = p)
  }
  y_int <- if (needs_int) {
    y_int_mat <- y_mat
    storage.mode(y_int_mat) <- "integer"
    y_int_mat
  } else {
    matrix(0L, nrow = n, ncol = p)
  }
  if (length(theta_anchor) == 1L) {
    theta_anchor_vec <- rep(as.double(theta_anchor), p)
  } else if (length(theta_anchor) == p) {
    theta_anchor_vec <- as.double(theta_anchor)
  } else {
    gdpar_abort(
      sprintf(
        paste(
          "'theta_anchor' must be a scalar (broadcast to length %d) or a",
          "vector of length %d; received length %d."
        ),
        p, p, length(theta_anchor)
      ),
      class = "gdpar_input_error",
      data = list(p = p, received_length = length(theta_anchor))
    )
  }
  X_mat <- if (d > 0L) design$X else matrix(0, nrow = n, ncol = 0L)
  use_dispersion_y <- as.integer(family$name == "gaussian")
  use_dispersion_phi <- as.integer(family$name == "neg_binomial_2")
  grp <- .resolve_group_id(group_id, n)
  list(
    n                  = n,
    p                  = p,
    family_id          = family$stan_id,
    use_a              = use_a,
    use_b              = use_b,
    use_W              = use_W,
    J_a_max            = as.integer(J_a_max),
    J_b_max            = as.integer(J_b_max),
    J_a_per_k          = as.integer(J_a_per_k),
    J_b_per_k          = as.integer(J_b_per_k),
    Z_a                = Z_a_padded,
    Z_b                = Z_b_padded,
    dim_W              = as.integer(dim_W),
    d                  = as.integer(d),
    W_per_k_dim        = as.integer(W_per_k_dim),
    X                  = X_mat,
    y_real             = y_real,
    y_int              = y_int,
    theta_anchor       = theta_anchor_vec,
    use_dispersion_y   = use_dispersion_y,
    use_dispersion_phi = use_dispersion_phi,
    W_type_id          = W_stan$W_type_id,
    W_n_knots_full     = W_stan$W_n_knots_full,
    W_knots_full       = W_stan$W_knots_full,
    W_degree           = W_stan$W_degree,
    use_groups         = grp$use_groups,
    J_groups           = grp$J_groups,
    group_id           = grp$group_id
  )
}

#' Internal: read the Path 1 multivariate Stan template and substitute placeholders
#'
#' Multivariate counterpart of \code{generate_stan_code}. Loads the
#' \code{amm_distrib_multi.stan} template shipped with the package and
#' replaces the prior and parametrization placeholders. The additive
#' scaling rule differs from the univariate path: the multivariate
#' codegen emits \code{" * sigma_a[k]"} in NCP (per-coordinate scale)
#' or \code{segment()}-based per-coordinate CP priors, controlled by
#' the logical vector \code{cp_a_per_k}. Three branches are wired
#' end-to-end via \code{\link{generate_a_blocks_multi}}: uniform NCP
#' (bit-exact with the H.1 codegen), uniform CP, and the mixed pattern
#' (caller injects \code{stan_data$cp_a_per_k_data} as an integer
#' vector, the template ramifies at runtime). The mixed pattern is the
#' Phase H.2 wiring of preflight CP/NCP per-\eqn{k} (Block 5.2).
#'
#' The modulating component placeholders \code{\{\{W_SCALE\}\}} and
#' \code{\{\{W_PRIOR\}\}} remain identical to the univariate path
#' because \code{W} is a globally shared parameter (single
#' \code{sigma_W[1]} regardless of \eqn{p}).
#'
#' Caching is delegated to \code{cmdstanr}, which hashes the source
#' string and reuses the compiled binary across calls when the source
#' is unchanged. CP and NCP variants therefore compile separately and
#' cache independently.
#'
#' @param prior An object of class \code{gdpar_prior}.
#' @param cp_a Logical scalar. Broadcast to a one-element pseudo-vector
#'   when \code{cp_a_per_k} is \code{NULL} (legacy mode). For
#'   per-coordinate control, pass \code{cp_a_per_k} explicitly.
#' @param cp_W Logical scalar. If \code{TRUE}, the modulating
#'   component \code{W} is sampled in the centered parametrization.
#'   If \code{FALSE} (default), non-centered.
#' @return A character scalar containing the Stan model source ready
#'   for \code{cmdstanr::cmdstan_model}.
#' @keywords internal
#' @noRd
generate_stan_code_multi <- function(prior, cp_a = FALSE, cp_W = FALSE,
                                     cp_a_per_k = NULL,
                                     template_name = "amm_distrib_multi.stan") {
  assert_inherits(prior, "gdpar_prior", "prior")
  # Translate legacy template names to their canonical-piece equivalents
  # in inst/stan/_canonical_pieces/. The FB default amm_distrib_multi.stan
  # was relocated in Sub-sub-fase 9.3.a (Sesion B9.3, 2026-05-27); the
  # EB-side multi templates (amm_eb_marginal_multi.stan +
  # amm_eb_conditional_multi.stan, canonized in Sub-fase 8.6.C decision
  # D34 to share the placeholder set bit-for-bit with this generator)
  # were relocated under the same canonical convention in Sub-bloque
  # 9.3.a colateral L.iv.A.2 (Bloque 9, Sesion B9.6, 2026-05-27): the
  # body pieces moved to _canonical_pieces/ and their inline helpers
  # block was reduced to the `// {{CANONICAL_HELPERS}}` placeholder
  # reusing the G.iv infrastructure of B9.4. Bit-exact Stan semantics
  # preservation verified token-wise.
  effective_template <- switch(template_name,
    "amm_distrib_multi.stan"      = "amm_canonical_pmulti.stan",
    "amm_eb_marginal_multi.stan"  = "amm_canonical_eb_marginal_multi.stan",
    "amm_eb_conditional_multi.stan" = "amm_canonical_eb_conditional_multi.stan",
    template_name
  )
  .gdpar_emit_canonical_stan(list(
    p_class = "pmulti",
    prior = prior,
    cp_a = cp_a,
    cp_W = cp_W,
    cp_a_per_k = cp_a_per_k,
    template_name = effective_template
  ))
}


#' Internal: generate the data/tp/model blocks for the a component
#'
#' Dispatch on the pattern of \code{cp_a_per_k}:
#' \itemize{
#'   \item All FALSE (NCP uniform): bit-exact with H.1 template.
#'     The transformed-parameters and model snippets are
#'     byte-identical to the H.1 substitution of
#'     \code{\{\{A_SCALE\}\} = " * sigma_a[k]"} and
#'     \code{\{\{A_PRIOR\}\} = "normal(0, 1)"}, and
#'     \code{\{\{DATA_CP_A_PER_K_DECL\}\}} is empty.
#'   \item All TRUE (CP uniform): a single Stan-side
#'     \code{for (k in 1:p)} loop with \code{segment()}-based
#'     per-coordinate CP priors; no extra data declaration.
#'   \item Mixed (some CP, some NCP): declares
#'     \code{array[p] int<lower=0, upper=1> cp_a_per_k_data} in
#'     the \code{data} block and uses a runtime ternary in
#'     transformed parameters and an \code{if/else} per coord in
#'     model. The caller must inject
#'     \code{stan_data$cp_a_per_k_data} as an integer vector.
#' }
#'
#' @param cp_a_per_k Logical vector. Length corresponds to \code{p}
#'   (the multivariate dimension).
#' @return Named list with character scalars \code{data_decl},
#'   \code{tp_block}, \code{model_block}.
#' @keywords internal
#' @noRd
generate_a_blocks_multi <- function(cp_a_per_k) {
  all_ncp <- !any(cp_a_per_k)
  all_cp  <- all(cp_a_per_k)
  mixed   <- !all_ncp && !all_cp
  if (all_ncp) {
    # Byte-identical (modulo placeholders) to the H.1 template
    # source after substitution of {{A_SCALE}} and {{A_PRIOR}}.
    return(list(
      data_decl   = "",
      tp_block    = paste(
        "  if (use_a == 1) {",
        "    for (k in 1:p) {",
        "      if (J_a_per_k[k] > 0) {",
        "        for (j in 1:J_a_free[k]) {",
        "          a_coef[k][j] = a_raw[a_raw_offset[k] + j] * sigma_a[k];",
        "        }",
        "      }",
        "    }",
        "  }",
        sep = "\n"
      ),
      model_block = "    a_raw ~ normal(0, 1);"
    ))
  }
  if (all_cp) {
    return(list(
      data_decl   = "",
      tp_block    = paste(
        "  if (use_a == 1) {",
        "    for (k in 1:p) {",
        "      if (J_a_per_k[k] > 0) {",
        "        for (j in 1:J_a_free[k]) {",
        "          a_coef[k][j] = a_raw[a_raw_offset[k] + j];",
        "        }",
        "      }",
        "    }",
        "  }",
        sep = "\n"
      ),
      model_block = paste(
        "    for (k in 1:p) {",
        "      if (J_a_free[k] > 0) {",
        "        segment(a_raw, a_raw_offset[k] + 1, J_a_free[k])",
        "          ~ normal(0, sigma_a[k]);",
        "      }",
        "    }",
        sep = "\n"
      )
    ))
  }
  # Mixed: per-k condition resolved at runtime via the data array
  # cp_a_per_k_data. Caller injects stan_data$cp_a_per_k_data.
  # Leading "\n" ensures the data declaration lands on its own line
  # (the placeholder is inline in the template to preserve bit-exact
  # paridad in the uniform-NCP branch, where data_decl is empty).
  list(
    data_decl   = "\n  array[p] int<lower=0, upper=1> cp_a_per_k_data;",
    tp_block    = paste(
      "  if (use_a == 1) {",
      "    for (k in 1:p) {",
      "      if (J_a_per_k[k] > 0) {",
      "        real scale_k = (cp_a_per_k_data[k] == 1) ? 1.0 : sigma_a[k];",
      "        for (j in 1:J_a_free[k]) {",
      "          a_coef[k][j] = a_raw[a_raw_offset[k] + j] * scale_k;",
      "        }",
      "      }",
      "    }",
      "  }",
      sep = "\n"
    ),
    model_block = paste(
      "    for (k in 1:p) {",
      "      if (J_a_free[k] > 0) {",
      "        if (cp_a_per_k_data[k] == 1) {",
      "          segment(a_raw, a_raw_offset[k] + 1, J_a_free[k])",
      "            ~ normal(0, sigma_a[k]);",
      "        } else {",
      "          segment(a_raw, a_raw_offset[k] + 1, J_a_free[k])",
      "            ~ normal(0, 1);",
      "        }",
      "      }",
      "    }",
      sep = "\n"
    )
  )
}

#' Internal: assemble the data block passed to the K-individual Stan template
#'
#' Multi-parametric counterpart of \code{assemble_stan_data} for the
#' regime \eqn{K > 1} with \eqn{p = 1}. Consumes the per-slot ragged
#' design returned by \code{.build_amm_design_K()} together with the
#' family that has had every K-individual slot promoted to
#' \code{scope = "per_observation"} by
#' \code{.gdpar_promote_scope_per_observation()} and produces the data
#' list expected by \code{amm_distrib_K.stan}.
#'
#' Sub-decision 1 of Unit 3 paso 3 (2026-05-20): the likelihood is a
#' distributional regression with canonical links per slot. Sub-phase
#' 8.3.4 closes Beta (\code{stan_id = 5}) and Gamma (\code{stan_id = 6})
#' alongside the Unit 3 baseline of Gaussian (\code{stan_id = 1}) and
#' neg_binomial_2 (\code{stan_id = 3}); sub-phase 8.3.5a closes
#' Student-t (\code{stan_id = 8}) as the K = 3 entry point. The
#' assembler delegates the family-id dispatch to the template via the
#' \code{family_id_k} array and asserts homogeneity (all slots share
#' the same \code{stan_id}) per gotcha 39 of handoff 29. Tweedie
#' (K >= 3) lands in sub-phase 8.3.5b; mixtures (ZIP / ZINB / Hurdle)
#' in 8.3.6.
#'
#' Sub-decision 2 of Unit 3 paso 3: the scope of W is global. The
#' assembler emits a single \code{use_W} flag and a single \code{dim_W}
#' derived from the first slot that declares W; \code{use_W = 1} when
#' any slot declares W. All slots that declare W must share the same
#' \code{W_basis} object (enforced by the canonical constructor
#' \code{.gdpar_formula_set_to_amm_spec_list()} which assigns the
#' external W argument to every slot whose RHS contains \code{W()}).
#'
#' Sub-decision 3 of Unit 3 paso 3: \code{theta_anchor_K} is per-slot.
#' This assembler expects \code{theta_anchor_K} pre-resolved into a
#' numeric vector of length K with names matching the slot names; the
#' caller (\code{.gdpar_K()}) handles broadcasting / validation.
#'
#' Gotcha 40 of handoff 29: \code{use_dispersion_y_k} and
#' \code{use_dispersion_phi_k} are emitted as integer vectors of length
#' K but are uniformly zero in Unit 3 because the dispersion slot
#' (\code{sigma} for Gaussian, \code{phi} for neg_binomial_2) is
#' K-individual and consumes its own AMM. The flags remain for future
#' extension to mixed (K-individual + population) family patterns in
#' sub-phase 8.3.4 onwards.
#'
#' @param design_K List returned by \code{.build_amm_design_K()}.
#' @param family A \code{gdpar_family} object whose \code{param_specs}
#'   have been promoted to \code{scope = "per_observation"} on the
#'   slots declared K-individual by the formula set.
#' @param amm_list_canonical Named list of \code{amm_spec} objects of
#'   length K. Used to extract the per-slot \code{use_a}, \code{use_b},
#'   \code{W} flags consistent with the design matrices.
#' @param y Numeric or integer vector of outcomes (length n).
#' @param theta_anchor_K Numeric vector of length K with the per-slot
#'   anchor values, names matching \code{design_K$slot_names}.
#' @param group_id Optional integer vector of length n carrying the
#'   group index of each observation. \code{NULL} (default) selects
#'   the single-group regime.
#'
#' @return Named list with the entries expected by
#'   \code{amm_distrib_K.stan}.
#' @keywords internal
#' @noRd
.assemble_stan_data_K <- function(design_K, family, amm_list_canonical, y,
                                   theta_anchor_K, group_id = NULL,
                                   family_id_k_vector = NULL) {
  if (!is.list(design_K) || is.null(design_K$Z_a_k_list)) {
    gdpar_abort(
      "Internal error: .assemble_stan_data_K expected a design list from .build_amm_design_K().",
      class = "gdpar_internal_error"
    )
  }
  K <- as.integer(design_K$K)
  if (K < 2L) {
    gdpar_abort(
      sprintf(
        "Internal error: .assemble_stan_data_K requires K >= 2; received K = %d.",
        K
      ),
      class = "gdpar_internal_error"
    )
  }
  if (length(amm_list_canonical) != K) {
    gdpar_abort(
      sprintf(
        "Internal error: amm_list_canonical has %d entries but design_K reports K = %d.",
        length(amm_list_canonical), K
      ),
      class = "gdpar_internal_error"
    )
  }
  if (inherits(family, "gdpar_family_multi")) {
    gdpar_abort(
      paste(
        "Internal error: the K-individual assembler does not accept",
        "gdpar_family_multi. The K > 1 + p > 1 regime is queued for Block 8.1."
      ),
      class = "gdpar_internal_error",
      data = list(received_class = class(family))
    )
  }
  assert_inherits(family, "gdpar_family", "family")
  if (is.null(family$stan_id) || is.na(family$stan_id)) {
    gdpar_abort(
      paste(
        "Internal error: family$stan_id is NA. Custom families with",
        "stan_id = NA cannot route through the built-in Stan template;",
        "Stan-side custom-family wiring (D-A3.B) lands in sub-phase 8.3.4."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(family_name = family$name)
    )
  }
  stan_id <- as.integer(family$stan_id)
  if (!stan_id %in% c(1L, 3L, 5L, 6L, 7L, 8L, 9L,
                      10L, 11L, 12L, 13L)) {
    gdpar_abort(
      sprintf(
        paste(
          "Family '%s' (stan_id = %d) is not wired in the K > 1 template",
          "amm_distrib_K.stan. Built-in K > 1 branches after sub-phase",
          "8.3.6: Gaussian (stan_id = 1), neg_binomial_2 (stan_id = 3),",
          "Beta (stan_id = 5), Gamma (stan_id = 6), the canonical",
          "custom-family pattern lognormal_loc_scale (stan_id = 7),",
          "Student-t (stan_id = 8), Tweedie (stan_id = 9), ZIP",
          "(stan_id = 10), ZINB (stan_id = 11), hurdle_poisson",
          "(stan_id = 12), and hurdle_neg_binomial_2 (stan_id = 13).",
          "Further custom-family patterns extend the registry returned",
          "by .gdpar_K_custom_patterns()."
        ),
        family$name, stan_id
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(family_name = family$name, stan_id = stan_id)
    )
  }
  n <- length(y)
  if (is.matrix(y) && ncol(y) > 1L) {
    gdpar_abort(
      paste(
        "When K > 1 with p = 1, the outcome 'y' must be a length-n vector",
        "(or a single-column matrix). Multivariate outcomes (p > 1) belong",
        "to the K = 1 multi-coord path; the combined K > 1 + p > 1 regime",
        "is queued for Block 8.1."
      ),
      class = "gdpar_input_error",
      data = list(received_ncol = ncol(y))
    )
  }
  y_vec <- if (is.matrix(y)) as.numeric(y[, 1L]) else as.numeric(y)
  if (length(y_vec) != n) {
    n <- length(y_vec)
  }
  grp <- .resolve_group_id(group_id, n)
  use_a_k <- vapply(amm_list_canonical,
                    function(a) as.integer(!is.null(a$a)),
                    integer(1L))
  use_b_k <- vapply(amm_list_canonical,
                    function(a) as.integer(!is.null(a$b)),
                    integer(1L))
  W_list <- lapply(amm_list_canonical, function(a) a$W)
  any_W <- any(!vapply(W_list, is.null, logical(1L)))
  use_W <- as.integer(any_W)
  J_a_per_k <- vapply(design_K$Z_a_k_list, ncol, integer(1L))
  J_b_per_k <- vapply(design_K$Z_b_k_list, ncol, integer(1L))
  J_a_max <- if (length(J_a_per_k) > 0L) max(J_a_per_k) else 0L
  J_b_max <- if (length(J_b_per_k) > 0L) max(J_b_per_k) else 0L
  pad_to <- function(z, target_cols, n_rows) {
    if (target_cols == 0L) {
      return(matrix(0, nrow = n_rows, ncol = 0L))
    }
    if (ncol(z) >= target_cols) {
      return(z)
    }
    cbind(z, matrix(0, nrow = nrow(z), ncol = target_cols - ncol(z)))
  }
  Z_a_padded <- lapply(design_K$Z_a_k_list, pad_to,
                       target_cols = J_a_max, n_rows = n)
  Z_b_padded <- lapply(design_K$Z_b_k_list, pad_to,
                       target_cols = J_b_max, n_rows = n)
  d <- ncol(design_K$X)
  if (use_W == 1L) {
    W_obj <- NULL
    for (k in seq_len(K)) {
      if (!is.null(W_list[[k]])) {
        W_obj <- W_list[[k]]
        break
      }
    }
    dim_W <- if (!is.null(W_obj$dim)) as.integer(W_obj$dim) else 0L
  } else {
    dim_W <- 0L
    W_obj <- NULL
  }
  W_stan <- .gdpar_resolve_W_stan_data(
    W = W_obj, use_W = use_W, y_vec = y_vec, family = family
  )
  needs_real <- stan_id %in% c(1L, 5L, 6L, 7L, 8L, 9L)
  needs_int <- stan_id %in% c(2L, 3L, 4L, 10L, 11L, 12L, 13L)
  if (stan_id == 5L) {
    if (any(!is.finite(y_vec)) || any(y_vec <= 0) || any(y_vec >= 1)) {
      gdpar_abort(
        paste(
          "Beta family requires the outcome y to lie strictly in (0, 1).",
          "Observed range falls on or outside the open unit interval;",
          "boundary observations are not handled by the canonical",
          "beta_proportion likelihood. Rescale or trim before fitting."
        ),
        class = "gdpar_input_error",
        data = list(family_name = family$name,
                    y_range = range(y_vec, na.rm = TRUE))
      )
    }
  }
  if (stan_id %in% c(6L, 7L)) {
    if (any(!is.finite(y_vec)) || any(y_vec <= 0)) {
      label <- if (stan_id == 6L) "Gamma" else "Lognormal_loc_scale"
      gdpar_abort(
        sprintf(
          paste(
            "%s family requires the outcome y to be strictly positive.",
            "Observed minimum is non-positive or non-finite; rescale or",
            "shift before fitting."
          ),
          label
        ),
        class = "gdpar_input_error",
        data = list(family_name = family$name,
                    y_range = range(y_vec, na.rm = TRUE))
      )
    }
  }
  if (stan_id == 9L) {
    if (any(!is.finite(y_vec)) || any(y_vec < 0)) {
      gdpar_abort(
        paste(
          "Tweedie family requires the outcome y to be non-negative.",
          "The canonical compound Poisson-gamma regime (1 < p < 2)",
          "admits a point mass at zero plus a continuous positive body;",
          "observations strictly less than zero or non-finite cannot be",
          "modeled. Rescale or shift before fitting."
        ),
        class = "gdpar_input_error",
        data = list(family_name = family$name,
                    y_range = range(y_vec, na.rm = TRUE))
      )
    }
  }
  if (stan_id %in% c(10L, 11L, 12L, 13L)) {
    label <- switch(as.character(stan_id),
      "10" = "ZIP",
      "11" = "ZINB",
      "12" = "hurdle_poisson",
      "13" = "hurdle_neg_binomial_2"
    )
    if (any(!is.finite(y_vec)) || any(y_vec < 0) ||
        any(y_vec != floor(y_vec))) {
      gdpar_abort(
        sprintf(
          paste(
            "%s family requires the outcome y to be a non-negative",
            "integer-valued count. Observations are either non-finite,",
            "strictly less than zero, or non-integer. Rescale, round or",
            "trim before fitting; mixture and hurdle models are count",
            "regression families, not continuous-with-mass."
          ),
          label
        ),
        class = "gdpar_input_error",
        data = list(family_name = family$name,
                    y_range = range(y_vec, na.rm = TRUE))
      )
    }
  }
  y_real <- if (needs_real) as.double(y_vec) else rep(0.0, n)
  y_int <- if (needs_int) as.integer(y_vec) else rep(0L, n)
  if (length(theta_anchor_K) != K) {
    gdpar_abort(
      sprintf(
        paste(
          "Internal error: theta_anchor_K must have length K = %d; received %d.",
          "Use the caller (.gdpar_K) to broadcast scalar / named-vector inputs."
        ),
        K, length(theta_anchor_K)
      ),
      class = "gdpar_internal_error",
      data = list(K = K, received_length = length(theta_anchor_K))
    )
  }
  X_mat <- if (d > 0L) design_K$X else matrix(0, nrow = n, ncol = 0L)
  if (is.null(family_id_k_vector)) {
    family_id_k <- rep(stan_id, K)
  } else {
    if (!is.numeric(family_id_k_vector) ||
        length(family_id_k_vector) != K ||
        any(is.na(family_id_k_vector))) {
      gdpar_abort(
        sprintf(
          paste(
            "Internal error: family_id_k_vector must be an integer vector",
            "of length K = %d with no NA; received length %d."
          ),
          K, length(family_id_k_vector)
        ),
        class = "gdpar_internal_error",
        data = list(K = K, received = family_id_k_vector)
      )
    }
    if (as.integer(family_id_k_vector[1L]) != stan_id) {
      gdpar_abort(
        sprintf(
          paste(
            "Internal error: heterogeneous family vector slot 1 stan_id (%d)",
            "does not match the location family stan_id (%d). The dispatcher",
            "must pass the location-slot family as `family` and the vector",
            "of per-slot stan_ids (including slot 1) as `family_id_k_vector`."
          ),
          as.integer(family_id_k_vector[1L]), stan_id
        ),
        class = "gdpar_internal_error",
        data = list(slot1_vector = family_id_k_vector[1L],
                    family_stan_id = stan_id)
      )
    }
    family_id_k <- as.integer(family_id_k_vector)
    allowed_ids <- c(1L, 3L, 5L, 6L, 7L, 8L, 9L,
                     10L, 11L, 12L, 13L)
    bad_slots <- which(!(family_id_k %in% allowed_ids))
    if (length(bad_slots) > 0L) {
      gdpar_abort(
        sprintf(
          paste(
            "family_id_k_vector contains unsupported stan_id(s) at slot(s)",
            "%s: %s. Built-in family stan_ids wired after sub-phase 8.3.7:",
            "1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13."
          ),
          paste(bad_slots, collapse = ", "),
          paste(family_id_k[bad_slots], collapse = ", ")
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(bad_slots = bad_slots,
                    bad_ids = family_id_k[bad_slots])
      )
    }
  }
  inv_link_id_per_slot <- .gdpar_compute_inv_link_id_per_slot(
    family_id_k, family
  )
  list(
    n                    = n,
    K                    = K,
    family_id_k          = family_id_k,
    inv_link_id_per_slot = inv_link_id_per_slot,
    use_a_k              = as.integer(use_a_k),
    use_b_k              = as.integer(use_b_k),
    use_W                = use_W,
    J_a_per_k            = as.integer(J_a_per_k),
    J_b_per_k            = as.integer(J_b_per_k),
    J_a_max              = as.integer(J_a_max),
    J_b_max              = as.integer(J_b_max),
    Z_a_k                = Z_a_padded,
    Z_b_k                = Z_b_padded,
    dim_W                = as.integer(dim_W),
    d                    = as.integer(d),
    X                    = X_mat,
    y_real               = y_real,
    y_int                = y_int,
    theta_anchor_K       = as.double(theta_anchor_K),
    use_dispersion_y_k   = as.integer(rep(0L, K)),
    use_dispersion_phi_k = as.integer(rep(0L, K)),
    W_type_id            = W_stan$W_type_id,
    W_n_knots_full       = W_stan$W_n_knots_full,
    W_knots_full         = W_stan$W_knots_full,
    W_degree             = W_stan$W_degree,
    use_groups           = grp$use_groups,
    J_groups             = grp$J_groups,
    group_id             = grp$group_id
  )
}

#' Internal: generate the data/tp/model blocks for the a component (K-individual)
#'
#' K-individual counterpart of \code{generate_a_blocks_multi}. The
#' structure is identical up to the index swap: the multivariate path
#' iterates over the AMM coordinate \eqn{k = 1, \dots, p}; the
#' K-individual path iterates over the family slot \eqn{k = 1, \dots, K}.
#' The three branches (NCP uniform / CP uniform / mixed) mirror the
#' multivariate path, including the flat-pack of \code{a_raw} via
#' per-slot \code{a_raw_offset}.
#'
#' @param cp_a_per_K Logical vector of length K. Length corresponds to
#'   the number of K-individual slots.
#' @return Named list with character scalars \code{data_decl},
#'   \code{tp_block}, \code{model_block}.
#' @keywords internal
#' @noRd
generate_a_blocks_K <- function(cp_a_per_K) {
  all_ncp <- !any(cp_a_per_K)
  all_cp  <- all(cp_a_per_K)
  mixed   <- !all_ncp && !all_cp
  if (all_ncp) {
    return(list(
      data_decl   = "",
      tp_block    = paste(
        "  if (any_use_a == 1) {",
        "    for (k in 1:K) {",
        "      if (use_a_k[k] == 1 && J_a_per_k[k] > 0) {",
        "        for (j in 1:J_a_free[k]) {",
        "          a_coef_k[k][j] = a_raw[a_raw_offset[k] + j] * sigma_a_k[sigma_a_idx[k]];",
        "        }",
        "      }",
        "    }",
        "  }",
        sep = "\n"
      ),
      model_block = "    a_raw ~ normal(0, 1);"
    ))
  }
  if (all_cp) {
    return(list(
      data_decl   = "",
      tp_block    = paste(
        "  if (any_use_a == 1) {",
        "    for (k in 1:K) {",
        "      if (use_a_k[k] == 1 && J_a_per_k[k] > 0) {",
        "        for (j in 1:J_a_free[k]) {",
        "          a_coef_k[k][j] = a_raw[a_raw_offset[k] + j];",
        "        }",
        "      }",
        "    }",
        "  }",
        sep = "\n"
      ),
      model_block = paste(
        "    for (k in 1:K) {",
        "      if (J_a_free[k] > 0) {",
        "        segment(a_raw, a_raw_offset[k] + 1, J_a_free[k])",
        "          ~ normal(0, sigma_a_k[sigma_a_idx[k]]);",
        "      }",
        "    }",
        sep = "\n"
      )
    ))
  }
  list(
    data_decl   = "\n  array[K] int<lower=0, upper=1> cp_a_per_K_data;",
    tp_block    = paste(
      "  if (any_use_a == 1) {",
      "    for (k in 1:K) {",
      "      if (use_a_k[k] == 1 && J_a_per_k[k] > 0) {",
      "        real scale_k = (cp_a_per_K_data[k] == 1) ? 1.0 : sigma_a_k[sigma_a_idx[k]];",
      "        for (j in 1:J_a_free[k]) {",
      "          a_coef_k[k][j] = a_raw[a_raw_offset[k] + j] * scale_k;",
      "        }",
      "      }",
      "    }",
      "  }",
      sep = "\n"
    ),
    model_block = paste(
      "    for (k in 1:K) {",
      "      if (J_a_free[k] > 0) {",
      "        if (cp_a_per_K_data[k] == 1) {",
      "          segment(a_raw, a_raw_offset[k] + 1, J_a_free[k])",
      "            ~ normal(0, sigma_a_k[sigma_a_idx[k]]);",
      "        } else {",
      "          segment(a_raw, a_raw_offset[k] + 1, J_a_free[k])",
      "            ~ normal(0, 1);",
      "        }",
      "      }",
      "    }",
      sep = "\n"
    )
  )
}

#' Internal: build the THETA_REF_PRIOR_BLOCK substitution per family
#'
#' Returns the Stan code block that sets the anchor priors per slot
#' for the K-individual distributional regression template. The
#' canonical (default) case keeps the vectorized prior used by every
#' family with homogeneous slot scales (Gaussian, neg_binomial_2,
#' Beta, Gamma, lognormal_loc_scale, Student-t). Tweedie (stan_id 9)
#' carries a structurally tight power slot \code{p in (1.01, 1.99)}
#' that cannot share the vectorized canonical prior; per decision E7
#' of D5 (sub-phase 8.3.5b, 2026-05-21), the block applies the
#' canonical prior to slots 1:2 via slice notation and the canonical
#' \code{uniform(1.01, 1.99)} prior to slot 3 separately.
#'
#' The substituted block still contains the \code{\{\{PRIOR_THETA_REF\}\}}
#' sub-placeholder which the calling \code{generate_stan_code_K}
#' resolves on the next pass; the helper therefore relies on the
#' caller ordering \code{THETA_REF_PRIOR_BLOCK} before
#' \code{PRIOR_THETA_REF} in its replacement list.
#'
#' @param family A \code{gdpar_family} object or \code{NULL}; the
#'   latter falls back to the default vectorized block.
#' @return A character scalar with the Stan code for the anchor-prior
#'   block (still containing the \code{PRIOR_THETA_REF} placeholder).
#' @keywords internal
#' @noRd
.gdpar_build_theta_ref_prior_block <- function(family) {
  stan_id <- if (is.null(family) || is.null(family$stan_id) ||
                 is.na(family$stan_id)) {
    NA_integer_
  } else {
    as.integer(family$stan_id)
  }
  default_block <- paste(
    "if (use_groups == 1) {",
    "    mu_theta_ref_k[1] ~ {{PRIOR_THETA_REF}};",
    "    sigma_theta_ref_k[1] ~ {{PRIOR_SIGMA_THETA_REF}};",
    "    for (g in 1:J_groups) {",
    "      theta_ref_k[g] ~ normal(mu_theta_ref_k[1], sigma_theta_ref_k[1]);",
    "    }",
    "  } else {",
    "    theta_ref_k[1] ~ {{PRIOR_THETA_REF}};",
    "  }",
    sep = "\n  "
  )
  if (is.na(stan_id) || stan_id != 9L) {
    return(default_block)
  }
  paste(
    "if (use_groups == 1) {",
    "    mu_theta_ref_k[1][1:2] ~ {{PRIOR_THETA_REF}};",
    "    mu_theta_ref_k[1][3] ~ uniform(1.01, 1.99);",
    "    sigma_theta_ref_k[1] ~ {{PRIOR_SIGMA_THETA_REF}};",
    "    for (g in 1:J_groups) {",
    "      theta_ref_k[g] ~ normal(mu_theta_ref_k[1], sigma_theta_ref_k[1]);",
    "    }",
    "  } else {",
    "    theta_ref_k[1][1:2] ~ {{PRIOR_THETA_REF}};",
    "    theta_ref_k[1][3] ~ uniform(1.01, 1.99);",
    "  }",
    sep = "\n  "
  )
}

#' Internal: build the THETA_REF_PRIOR_BLOCK substitution per family
#' for the multi-parametric multivariate (KxP) regime
#'
#' Parallel to \code{.gdpar_build_theta_ref_prior_block} but emits the
#' anchor prior block for the FB Path C template
#' \code{amm_canonical_pmulti_KxP.stan} (Sub-bloque 9.3.d, Sesion B9.5,
#' 2026-05-27). The anchor container is
#' \code{array[J_groups, K] vector[p] theta_ref_kp} instead of the
#' \code{array[J_groups] vector[K] theta_ref_k} of the K piece; the
#' hierarchical hyperprior is per-slot-per-coord
#' \code{(mu_theta_ref_kp[g, k], sigma_theta_ref_kp[g, k])}. The block
#' iterates over the K slots so the substitution does not have to
#' explicitly know K (which is data, not a literal).
#'
#' Tweedie (stan_id 9) carries the structurally tight slot p in
#' (1.01, 1.99); per E7 of D5 (sub-phase 8.3.5b) the block applies the
#' canonical vectorized prior to slots 1:2 of every coord and the
#' \code{uniform(1.01, 1.99)} prior to slot 3 of every coord
#' separately.
#'
#' @param family A \code{gdpar_family} object or \code{NULL}; the
#'   latter falls back to the default vectorized block.
#' @return A character scalar with the Stan code for the anchor-prior
#'   block (still containing the \code{PRIOR_THETA_REF} and
#'   \code{PRIOR_SIGMA_THETA_REF} placeholders for the next pass).
#' @keywords internal
#' @noRd
.gdpar_build_theta_ref_prior_block_KxP <- function(family) {
  stan_id <- if (is.null(family) || is.null(family$stan_id) ||
                 is.na(family$stan_id)) {
    NA_integer_
  } else {
    as.integer(family$stan_id)
  }
  default_block <- paste(
    "if (use_groups == 1) {",
    "    for (k in 1:K) {",
    "      mu_theta_ref_kp[1, k] ~ {{PRIOR_THETA_REF}};",
    "      sigma_theta_ref_kp[1, k] ~ {{PRIOR_SIGMA_THETA_REF}};",
    "      for (g in 1:J_groups) {",
    "        theta_ref_kp[g, k] ~ normal(mu_theta_ref_kp[1, k],",
    "                                    sigma_theta_ref_kp[1, k]);",
    "      }",
    "    }",
    "  } else {",
    "    for (k in 1:K) {",
    "      theta_ref_kp[1, k] ~ {{PRIOR_THETA_REF}};",
    "    }",
    "  }",
    sep = "\n  "
  )
  if (is.na(stan_id) || stan_id != 9L) {
    return(default_block)
  }
  paste(
    "if (use_groups == 1) {",
    "    for (k in 1:K) {",
    "      if (k <= 2) {",
    "        mu_theta_ref_kp[1, k] ~ {{PRIOR_THETA_REF}};",
    "      } else {",
    "        for (jj in 1:p) {",
    "          mu_theta_ref_kp[1, k][jj] ~ uniform(1.01, 1.99);",
    "        }",
    "      }",
    "      sigma_theta_ref_kp[1, k] ~ {{PRIOR_SIGMA_THETA_REF}};",
    "      for (g in 1:J_groups) {",
    "        theta_ref_kp[g, k] ~ normal(mu_theta_ref_kp[1, k],",
    "                                    sigma_theta_ref_kp[1, k]);",
    "      }",
    "    }",
    "  } else {",
    "    for (k in 1:K) {",
    "      if (k <= 2) {",
    "        theta_ref_kp[1, k] ~ {{PRIOR_THETA_REF}};",
    "      } else {",
    "        for (jj in 1:p) {",
    "          theta_ref_kp[1, k][jj] ~ uniform(1.01, 1.99);",
    "        }",
    "      }",
    "    }",
    "  }",
    sep = "\n  "
  )
}

#' Internal: read the Path C (KxP) FB Stan template and substitute placeholders
#'
#' Thin wrapper around \code{.gdpar_emit_canonical_stan} for the
#' multi-parametric multivariate Full-Bayes regime (K >= 2 slots
#' crossed with p >= 2 coords), introduced in Sub-bloque 9.3.d
#' (Bloque 9, Sesion B9.5, 2026-05-27) under decisions I.iv lateral +
#' J.iv.A piece arquitectonica. The wrapper preserves the canonical
#' source-of-truth pattern of B9.3 (dispatcher reads piece, applies
#' substitutions) and is intentionally internal: the user-facing
#' guard \code{.gdpar_guard_multiparam_multivariate} still aborts
#' \code{gdpar()} on K > 1 + p > 1, deferring the public surface lift
#' to a future sub-phase. The wrapper exists so the bootstrap scripts
#' in \code{data-raw/bootstrap_fb_goldens_KxP.R} can emit and compile
#' the Path C template deterministically per the reserved canonical
#' seeds 91001..91004 documented in
#' \code{~/claude-setup/gdpar-block-9/DESIGN_9_3_D_PATH_C.md} §4.
#'
#' @param prior An object of class \code{gdpar_prior}.
#' @param cp_W Logical scalar. If \code{TRUE}, the modulating
#'   component W is sampled in the centered parametrization; if
#'   \code{FALSE} (default), non-centered.
#' @param family Optional \code{gdpar_family} object; when supplied
#'   drives the family-specific anchor-prior block emission via
#'   \code{.gdpar_build_theta_ref_prior_block_KxP()}. \code{NULL}
#'   (default) keeps the canonical vectorized block.
#' @return A character scalar containing the Stan model source ready
#'   for \code{cmdstanr::cmdstan_model}.
#' @keywords internal
#' @noRd
generate_stan_code_KxP_FB <- function(prior, cp_W = FALSE,
                                       family = NULL) {
  assert_inherits(prior, "gdpar_prior", "prior")
  .gdpar_emit_canonical_stan(list(
    p_class = "pmulti_KxP",
    prior   = prior,
    cp_W    = cp_W,
    family  = family
  ))
}

#' Internal: FB Path C fit harness for the K > 1 + p > 1 regime
#'
#' Sub-bloque 9.3.d colateral D69 (Bloque 9, Sesion B9.6, 2026-05-27)
#' under decision N.a canonized at session open. Closes the operational
#' side of the Path C bit-exact comparator template canonized
#' architecturally in B9.5 (decision J.iv.A piece dedicada
#' \code{amm_canonical_pmulti_KxP.stan} + thin wrapper
#' \code{generate_stan_code_KxP_FB} + helper
#' \code{.gdpar_build_theta_ref_prior_block_KxP}). Composes:
#' \enumerate{
#'   \item the canonical KxP design via \code{.build_amm_design_KxP()},
#'   \item the unified KxP assembler with \code{path = "FB"} (decision
#'     M.iv unified canonized B9.6), and
#'   \item \code{cmdstanr::cmdstan_model()$sample()} on the substituted
#'     KxP Stan source.
#' }
#'
#' The harness intentionally bypasses the user-facing
#' \code{.gdpar_guard_multiparam_multivariate} guard (which still
#' aborts \code{gdpar()} for \code{K > 1 + p > 1}) since the canonical
#' piece + dispatcher already exist; the public lift of the guard is
#' deferred to a future sub-fase post-Bloque 9.
#'
#' The function is internal to the package and is consumed exclusively
#' by the bootstrap script \code{data-raw/bootstrap_fb_goldens_KxP.R}
#' and the Tier 1 compare-path tests under
#' \code{tests/testthat/test-fb_goldens_9_3_d_KxP.R}. It is not exported
#' through \code{NAMESPACE}.
#'
#' @param y_matrix Numeric matrix \code{n} by \code{p} with the
#'   outcomes.
#' @param family A \code{gdpar_family} object whose \code{stan_id} is in
#'   the Path B set \code{{1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13}}.
#' @param amm_list Named list of length \code{K >= 2} of
#'   \code{amm_spec} objects each with \code{p >= 2} per-slot. Slots
#'   must share \code{p}.
#' @param data Data frame providing the covariates referenced by the
#'   per-slot \code{a} and \code{b} formulas of every \code{amm_spec}.
#' @param formula_rhs One-sided formula whose right-hand side enumerates
#'   the union of covariates referenced by the slot formulas plus the
#'   covariate(s) consumed by \code{W()} when present.
#' @param theta_anchor_kp Numeric matrix \code{K} by \code{p} with the
#'   per-slot per-coord anchors on the linear-predictor scale; defaults
#'   to a zero matrix.
#' @param iter_warmup,iter_sampling,chains,seed,refresh,... Forwarded
#'   to \code{cmdstanr::cmdstan_model()$sample()}.
#' @return A named list with the cmdstanr fit (\code{$fit}), the
#'   assembled Stan data (\code{$stan_data}), the family metadata, the
#'   sampling configuration and the canonical seed, suitable for
#'   bit-exact RDS production under
#'   \code{tests/testthat/data/golden_fb_KxP_*.rds}.
#' @keywords internal
#' @noRd
.gdpar_fb_KxP_fit <- function(y_matrix, family, amm_list, data,
                               formula_rhs,
                               theta_anchor_kp = NULL,
                               iter_warmup = 200L,
                               iter_sampling = 200L,
                               chains = 2L,
                               seed = 91001L,
                               refresh = 0L, ...) {
  assert_inherits(family, "gdpar_family", "family")
  if (!is.list(amm_list) || length(amm_list) < 2L) {
    gdpar_abort(
      "FB KxP fit requires amm_list to be a named list of K >= 2 amm_spec objects.",
      class = "gdpar_input_error",
      data = list(K = length(amm_list))
    )
  }
  K <- length(amm_list)
  p <- as.integer(amm_list[[1L]]$p %||% 1L)
  if (p < 2L) {
    gdpar_abort(
      "FB KxP fit requires p >= 2 per slot.",
      class = "gdpar_input_error",
      data = list(p = p)
    )
  }
  # Defensive: materialize W_basis if any slot carries one with unmaterialized dim.
  for (k in seq_along(amm_list)) {
    if (!is.null(amm_list[[k]]$W) &&
        (is.null(amm_list[[k]]$W$dim) ||
         is.na(amm_list[[k]]$W$dim))) {
      amm_list[[k]]$W <- materialize_W_basis(amm_list[[k]]$W, p = p)
    }
  }
  design_KxP <- .build_amm_design_KxP(amm_list, data, formula_rhs)
  if (is.null(theta_anchor_kp)) {
    theta_anchor_kp <- matrix(0, nrow = K, ncol = p)
  }
  if (!is.matrix(y_matrix) || nrow(y_matrix) != nrow(data) ||
      ncol(y_matrix) != p) {
    gdpar_abort(
      sprintf(paste(
        "y_matrix must be a numeric matrix of dimensions nrow(data) = %d",
        "by p = %d."
      ), nrow(data), p),
      class = "gdpar_input_error",
      data = list(nrow_y = if (is.matrix(y_matrix)) nrow(y_matrix) else NA,
                  ncol_y = if (is.matrix(y_matrix)) ncol(y_matrix) else NA)
    )
  }
  stan_data <- .assemble_stan_data_KxP(
    design_KxP, family, amm_list, y_matrix, theta_anchor_kp,
    path = "FB"
  )
  prior <- gdpar_prior()
  stan_src <- generate_stan_code_KxP_FB(prior, cp_W = FALSE, family = family)
  stan_path <- write_stan_to_tempfile(stan_src)
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    gdpar_abort(
      "Package 'cmdstanr' is required for .gdpar_fb_KxP_fit().",
      class = "gdpar_unsupported_feature_error"
    )
  }
  model <- cmdstanr::cmdstan_model(stan_path, compile = TRUE, quiet = TRUE)
  # Conservative cmdstanr init: the canonical KxP piece composes
  # Path A multivariate with Path B distributional regression; the
  # default cmdstanr init = 2 (unconstrained uniform [-2, 2]) maps to
  # inv_logit / exp blow-ups for the Path B families (beta, gamma,
  # student_t, ...) when their location / scale eta lands near the
  # boundary of the link support. Reuse the EB-side canonical helper
  # `.gdpar_eb_make_random_init_KxP()` (canonized in Sub-fase 8.6.D
  # Session 13c, 2026-05-26) which builds a chain-deterministic init
  # list with theta_ref_kp = small jitter, sigma_a_k = 0.1 + jitter,
  # sigma_b_k = 0.1 + jitter, a_raw = small jitter, c_b_kp_raw = small
  # jitter. Parameters not consumed by the FB-side extension
  # (sigma_W, W_raw when use_W = 0; sigma_y_pop_k, phi_pop_k when
  # use_dispersion_*_k all zero) are size-zero by construction so the
  # missing init entries are harmless.
  fb_init_fn <- function(chain_id = 1L) {
    .gdpar_eb_make_random_init_KxP(stan_data,
                                    seed_offset = as.integer(chain_id),
                                    base_seed = seed)
  }
  fit <- model$sample(
    data          = stan_data,
    iter_warmup   = iter_warmup,
    iter_sampling = iter_sampling,
    chains        = chains,
    refresh       = refresh,
    seed          = seed,
    init          = fb_init_fn,
    ...
  )
  list(
    fit              = fit,
    stan_data        = stan_data,
    family_name      = family$name,
    family_stan_id   = as.integer(family$stan_id),
    K                = as.integer(K),
    p                = as.integer(p),
    n                = as.integer(nrow(data)),
    seed             = as.integer(seed),
    iter_warmup      = as.integer(iter_warmup),
    iter_sampling    = as.integer(iter_sampling),
    chains           = as.integer(chains),
    cmdstan_version  = .gdpar_cmdstan_version()
  )
}

#' Internal: read the K-individual Stan template and substitute placeholders
#'
#' K-individual counterpart of \code{generate_stan_code_multi}. Loads
#' the \code{amm_distrib_K.stan} template shipped with the package
#' (Sub-decision tecnica a.i, two-file materialization, 2026-05-20)
#' and replaces the prior, parametrization and a-block placeholders.
#' The W placeholders (\code{\{\{W_SCALE\}\}}, \code{\{\{W_PRIOR\}\}})
#' remain identical to the multivariate path because W is a globally
#' shared parameter (single \code{sigma_W[1]}) and the placeholders
#' substitute the same expressions across templates.
#'
#' Sub-phase 8.3.5b adds the \code{family} argument so the function
#' can choose the per-family \code{THETA_REF_PRIOR_BLOCK} substitution
#' (see \code{.gdpar_build_theta_ref_prior_block()}). A \code{NULL}
#' value preserves the legacy vectorized block for backward
#' compatibility with existing test call sites.
#'
#' @param prior An object of class \code{gdpar_prior}.
#' @param cp_a Logical scalar. Broadcast to a one-element pseudo-vector
#'   when \code{cp_a_per_K} is \code{NULL} (legacy mode); per-slot
#'   control via \code{cp_a_per_K}.
#' @param cp_W Logical scalar. If \code{TRUE}, the modulating component
#'   W is sampled in the centered parametrization; if \code{FALSE}
#'   (default), non-centered.
#' @param cp_a_per_K Logical vector of length K. \code{NULL} (default)
#'   falls back to \code{cp_a} broadcast.
#' @param family Optional \code{gdpar_family} object; when supplied
#'   drives the family-specific anchor-prior block emission. \code{NULL}
#'   (default) keeps the canonical vectorized block.
#' @return A character scalar containing the Stan model source ready
#'   for \code{cmdstanr::cmdstan_model}.
#' @keywords internal
#' @noRd
generate_stan_code_K <- function(prior, cp_a = FALSE, cp_W = FALSE,
                                  cp_a_per_K = NULL, family = NULL,
                                  template_name = "amm_distrib_K.stan") {
  assert_inherits(prior, "gdpar_prior", "prior")
  # Sub-bloque 9.3.a colateral L.iv.A (Bloque 9, Sesion B9.5,
  # 2026-05-27) + L.iv.A.2 (Sesion B9.6, 2026-05-27): translate legacy
  # template names to their canonical-piece equivalents in
  # inst/stan/_canonical_pieces/. The FB default amm_distrib_K.stan was
  # relocated in B9.5; the EB-side K templates
  # (amm_eb_marginal_K.stan + amm_eb_conditional_K.stan, canonized in
  # D34 of Sub-fase 8.6.C) were relocated under the same convention in
  # B9.6 cascade L heredada. The inline bspline_basis_eval +
  # apply_W_basis_diff helpers were reduced to the
  # `// {{CANONICAL_HELPERS}}` placeholder (G.iv pattern); the inline
  # Tweedie functions block + apply_inv_link_by_id stay K-specific and
  # remain inline. Bit-exact Stan semantics preservation verified
  # token-wise.
  effective_template <- switch(template_name,
    "amm_distrib_K.stan"           = "amm_canonical_distrib_K.stan",
    "amm_eb_marginal_K.stan"       = "amm_canonical_eb_marginal_K.stan",
    "amm_eb_conditional_K.stan"    = "amm_canonical_eb_conditional_K.stan",
    template_name
  )
  .gdpar_emit_canonical_stan(list(
    p_class       = "distrib_K",
    prior         = prior,
    cp_a          = cp_a,
    cp_W          = cp_W,
    cp_a_per_K    = cp_a_per_K,
    family        = family,
    template_name = effective_template
  ))
}

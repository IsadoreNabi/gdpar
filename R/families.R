#' Internal: build a parameter specification
#'
#' Constructs a \code{gdpar_param_spec} describing one structural
#' parameter of a family under the canonical contract decided in the
#' multi-parametric extension scoping (Block 8 Session 1, decision 1C):
#' the family is the list of marginal parameter specifications. Each
#' spec carries the parameter name, link, identifiability status under
#' Lemma 1B, support, canonical prior kind (decision 2b-iii) and scope
#' (decision 2a-iii).
#'
#' @keywords internal
#' @noRd
.gdpar_param_spec <- function(name, link, did_status,
                              did_condition, did_reference,
                              support, prior_canonical_kind,
                              scope, family_role) {
  link_funcs <- .gdpar_link_funcs(link)
  obj <- list(
    name                 = name,
    link                 = link,
    linkfun              = link_funcs$linkfun,
    inv_link             = link_funcs$inv_link,
    did_status           = did_status,
    did_condition        = did_condition,
    did_reference        = did_reference,
    support              = support,
    prior_canonical_kind = prior_canonical_kind,
    scope                = scope,
    family_role          = family_role
  )
  class(obj) <- c("gdpar_param_spec", "list")
  obj
}

#' Internal: known family_role values
#'
#' Returns the character vector of all family roles registered for
#' \code{param_spec} entries. Block 8 Session 1 decision 1C
#' (family is a list of param_specs) plus the Block 8.3.2 extension to
#' arbitrary K and structural role per parameter. Roles:
#' \itemize{
#'   \item \code{location}: position parameter (mu in any GLMM-like
#'     family; the per-observation primary parameter).
#'   \item \code{scale}: dispersion / scale parameter (sigma for
#'     Gaussian, phi for NB/Beta, phi for Tweedie).
#'   \item \code{shape}: shape parameter (Gamma shape, when relevant).
#'   \item \code{df}: degrees of freedom (Student-t).
#'   \item \code{mixture_pi}: mixture probability (ZIP / ZINB /
#'     Hurdle).
#'   \item \code{power}: power parameter (Tweedie p in (1, 2)).
#' }
#'
#' Used by the validator inside \code{.gdpar_param_spec()} and by the
#' codegen branches that emit family-specific Stan likelihood blocks.
#' @return Character vector.
#' @keywords internal
#' @noRd
.gdpar_known_family_roles <- function() {
  c("location", "scale", "shape", "df", "mixture_pi", "power")
}

#' Internal: link function factory
#'
#' Returns a list with \code{linkfun} and \code{inv_link} for a given
#' link name. Centralizes the link switch used by the family and
#' parameter spec constructors.
#'
#' @keywords internal
#' @noRd
.gdpar_link_funcs <- function(link) {
  inv_link <- switch(link,
    identity = function(eta) eta,
    log      = function(eta) exp(eta),
    logit    = function(eta) 1 / (1 + exp(-eta)),
    gdpar_abort(
      sprintf("Internal error: unsupported link '%s'.", link),
      class = "gdpar_internal_error"
    )
  )
  linkfun <- switch(link,
    identity = function(mu) mu,
    log      = function(mu) log(mu),
    logit    = function(mu) log(mu / (1 - mu)),
    gdpar_abort(
      sprintf("Internal error: unsupported link '%s'.", link),
      class = "gdpar_internal_error"
    )
  )
  list(linkfun = linkfun, inv_link = inv_link)
}

#' Internal: canonical param_spec list per built-in family
#'
#' Returns the list of \code{gdpar_param_spec} objects describing every
#' structural parameter the built-in family admits as \emph{eligible}
#' for an individual specification. The first entry carries the
#' per-observation location parameter; subsequent entries are
#' population-level auxiliaries in their canonical scope. K-individual
#' membership (which auxiliaries are promoted to per-observation) is
#' declared at the entry of \code{\link{gdpar}} via a
#' \code{\link{gdpar_formula_set}} or a named list of
#' \code{\link{amm_spec}}; the helper
#' \code{.gdpar_promote_scope_per_observation()} consumes that
#' declaration to return a family copy with promoted \code{scope}
#' fields, which the codegen consumes under the hybrid shape of
#' decision 2a-iii. The family is the registry of eligibles, not the
#' subset of individuals.
#'
#' @keywords internal
#' @noRd
.gdpar_family_param_specs_for <- function(name, link) {
  switch(name,
    gaussian = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "real_line",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "sigma",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "variance identifiable from data",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_sigma",
        scope                = "population",
        family_role          = "scale"
      )
    ),
    poisson = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      )
    ),
    neg_binomial_2 = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "phi",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical overdispersion sufficient to identify phi",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_phi",
        scope                = "population",
        family_role          = "scale"
      )
    ),
    bernoulli = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "unit_interval",
        prior_canonical_kind = "logit_p",
        scope                = "per_observation",
        family_role          = "location"
      )
    ),
    beta = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "unit_interval",
        prior_canonical_kind = "logit_p",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "phi",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical dispersion sufficient to identify phi",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_phi",
        scope                = "population",
        family_role          = "scale"
      )
    ),
    gamma = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "shape",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical dispersion sufficient to identify shape",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_shape",
        scope                = "population",
        family_role          = "shape"
      )
    ),
    student_t = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "real_line",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "sigma",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical residual dispersion sufficient to identify sigma",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_sigma",
        scope                = "population",
        family_role          = "scale"
      ),
      .gdpar_param_spec(
        name                 = "nu",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical tail thickness sufficient to identify degrees of freedom",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_nu",
        scope                = "population",
        family_role          = "df"
      )
    ),
    tweedie = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "phi",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical dispersion sufficient to identify phi",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "positive_real",
        prior_canonical_kind = "log_phi",
        scope                = "population",
        family_role          = "scale"
      ),
      .gdpar_param_spec(
        name                 = "p",
        link                 = "identity",
        did_status           = "holds_under_condition",
        did_condition        = "empirical balance between zero mass and continuous body sufficient to identify the power parameter",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B)",
        support              = "bounded_open",
        prior_canonical_kind = "power_p",
        scope                = "population",
        family_role          = "power"
      )
    ),
    zip = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Lambert (1992)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "pi",
        link                 = "logit",
        did_status           = "holds_under_condition",
        did_condition        = "empirical proportion of structural zeros distinct from sampling zeros of the count component sufficient to identify the zero-inflation probability",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Lambert (1992)",
        support              = "unit_interval",
        prior_canonical_kind = "logit_pi",
        scope                = "population",
        family_role          = "mixture_pi"
      )
    ),
    zinb = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Greene (1994)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "phi",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical overdispersion sufficient to identify phi",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Greene (1994)",
        support              = "positive_real",
        prior_canonical_kind = "log_phi",
        scope                = "population",
        family_role          = "scale"
      ),
      .gdpar_param_spec(
        name                 = "pi",
        link                 = "logit",
        did_status           = "holds_under_condition",
        did_condition        = "empirical proportion of structural zeros distinct from sampling zeros of the count component sufficient to identify the zero-inflation probability",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Greene (1994)",
        support              = "unit_interval",
        prior_canonical_kind = "logit_pi",
        scope                = "population",
        family_role          = "mixture_pi"
      )
    ),
    hurdle_poisson = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Mullahy (1986)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "pi",
        link                 = "logit",
        did_status           = "holds_under_condition",
        did_condition        = "empirical proportion of zeros sufficient to identify the hurdle probability",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Mullahy (1986)",
        support              = "unit_interval",
        prior_canonical_kind = "logit_pi",
        scope                = "population",
        family_role          = "mixture_pi"
      )
    ),
    hurdle_neg_binomial_2 = list(
      .gdpar_param_spec(
        name                 = "mu",
        link                 = link,
        did_status           = "holds",
        did_condition        = NA_character_,
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Mullahy (1986)",
        support              = "positive_real",
        prior_canonical_kind = "mu",
        scope                = "per_observation",
        family_role          = "location"
      ),
      .gdpar_param_spec(
        name                 = "phi",
        link                 = "log",
        did_status           = "holds_under_condition",
        did_condition        = "empirical overdispersion sufficient to identify phi",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Mullahy (1986)",
        support              = "positive_real",
        prior_canonical_kind = "log_phi",
        scope                = "population",
        family_role          = "scale"
      ),
      .gdpar_param_spec(
        name                 = "pi",
        link                 = "logit",
        did_status           = "holds_under_condition",
        did_condition        = "empirical proportion of zeros sufficient to identify the hurdle probability",
        did_reference        = "Block 1, Section 6.4 (Lemma 1B); Mullahy (1986)",
        support              = "unit_interval",
        prior_canonical_kind = "logit_pi",
        scope                = "population",
        family_role          = "mixture_pi"
      )
    ),
    gdpar_abort(
      sprintf(
        "Internal error: no canonical param_specs for family '%s'.",
        name
      ),
      class = "gdpar_internal_error"
    )
  )
}

#' Internal: apply a user-supplied D-ID override to a list of param_specs
#'
#' Implements the \code{did_override} plasticity argument of
#' \code{\link{gdpar_family}} (D2 of sub-phase 8.3.5a, 2026-05-21):
#' the user keeps the canonical likelihood and links of a built-in
#' family but replaces the identifiability descriptors of one or more
#' slots with their own declaration. Returns a new list of param_specs
#' with the overrides applied; the input is not mutated.
#'
#' Validation: \code{did_override} must be either \code{NULL} (no-op)
#' or a named list whose names are a subset of the slot names of the
#' family's param_specs. Each entry must be a list whose names are a
#' subset of \code{c("did_status", "did_condition", "did_reference")};
#' unknown fields raise an error. \code{did_status} must be one of
#' \code{c("holds", "holds_under_condition", "user_responsible")}.
#'
#' Why descriptor-level override and not surface duplication via the
#' custom_K registry: D2 of the sub-phase 8.3.5 scoping (2026-05-21)
#' rejected registering a built-in inside
#' \code{.gdpar_K_custom_patterns()} because that blurs the contract
#' boundary between "package canonizes the likelihood" and "user
#' declares D-ID". The override keeps the likelihood + links
#' first-class while letting the user adjust the identifiability
#' semantics for their design.
#'
#' @param param_specs List of \code{gdpar_param_spec} objects.
#' @param did_override Either \code{NULL} or a named list (see Details).
#' @param family_name Character scalar for error messages.
#' @return List of param_specs with the overrides applied.
#' @keywords internal
#' @noRd
.gdpar_apply_did_override <- function(param_specs, did_override,
                                       family_name) {
  if (is.null(did_override)) {
    return(param_specs)
  }
  if (!is.list(did_override) || is.null(names(did_override)) ||
      any(names(did_override) == "")) {
    gdpar_abort(
      "Argument 'did_override' must be NULL or a named list keyed by slot name.",
      class = "gdpar_input_error"
    )
  }
  slot_names <- vapply(param_specs, function(s) s$name, character(1L))
  unknown <- setdiff(names(did_override), slot_names)
  if (length(unknown) > 0L) {
    gdpar_abort(
      sprintf(
        paste(
          "did_override contains slot name(s) %s not present in family '%s'.",
          "Eligible slot names: %s."
        ),
        paste(sQuote(unknown), collapse = ", "),
        family_name,
        paste(sQuote(slot_names), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(unknown = unknown, eligible = slot_names,
                  family = family_name)
    )
  }
  allowed_fields   <- c("did_status", "did_condition", "did_reference")
  allowed_statuses <- c("holds", "holds_under_condition", "user_responsible")
  for (slot in names(did_override)) {
    entry <- did_override[[slot]]
    if (!is.list(entry) || is.null(names(entry)) ||
        any(names(entry) == "")) {
      gdpar_abort(
        sprintf(
          "did_override[['%s']] must be a named list with optional fields %s.",
          slot, paste(sQuote(allowed_fields), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(slot = slot)
      )
    }
    bad_fields <- setdiff(names(entry), allowed_fields)
    if (length(bad_fields) > 0L) {
      gdpar_abort(
        sprintf(
          "did_override[['%s']] contains unknown field(s): %s. Allowed: %s.",
          slot, paste(sQuote(bad_fields), collapse = ", "),
          paste(sQuote(allowed_fields), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(slot = slot, unknown = bad_fields,
                    allowed = allowed_fields)
      )
    }
    if (!is.null(entry$did_status)) {
      if (!is.character(entry$did_status) ||
          length(entry$did_status) != 1L ||
          !entry$did_status %in% allowed_statuses) {
        gdpar_abort(
          sprintf(
            "did_override[['%s']]$did_status must be one of %s.",
            slot,
            paste(sQuote(allowed_statuses), collapse = ", ")
          ),
          class = "gdpar_input_error",
          data = list(slot = slot, received = entry$did_status,
                      allowed = allowed_statuses)
        )
      }
    }
    if (!is.null(entry$did_condition)) {
      if (!is.character(entry$did_condition) ||
          length(entry$did_condition) != 1L) {
        gdpar_abort(
          sprintf(
            "did_override[['%s']]$did_condition must be a character scalar (NA_character_ for unconditional).",
            slot
          ),
          class = "gdpar_input_error",
          data = list(slot = slot)
        )
      }
    }
    if (!is.null(entry$did_reference)) {
      if (!is.character(entry$did_reference) ||
          length(entry$did_reference) != 1L) {
        gdpar_abort(
          sprintf(
            "did_override[['%s']]$did_reference must be a character scalar.",
            slot
          ),
          class = "gdpar_input_error",
          data = list(slot = slot)
        )
      }
    }
  }
  for (i in seq_along(param_specs)) {
    slot <- param_specs[[i]]$name
    entry <- did_override[[slot]]
    if (is.null(entry)) next
    if (!is.null(entry$did_status)) {
      param_specs[[i]]$did_status <- entry$did_status
    }
    if (!is.null(entry$did_condition)) {
      param_specs[[i]]$did_condition <- entry$did_condition
    }
    if (!is.null(entry$did_reference)) {
      param_specs[[i]]$did_reference <- entry$did_reference
    }
  }
  param_specs
}

#' Print method for gdpar_param_spec objects
#'
#' @param x An object of class \code{gdpar_param_spec}.
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return Invisibly returns \code{x}.
#' @keywords internal
#' @exportS3Method base::print
print.gdpar_param_spec <- function(x, ...) {
  cat("<gdpar_param_spec>\n")
  cat("  name                 : ", x$name, "\n", sep = "")
  cat("  link                 : ", x$link, "\n", sep = "")
  cat("  family_role          : ", x$family_role, "\n", sep = "")
  cat("  scope                : ", x$scope, "\n", sep = "")
  cat("  support              : ", x$support, "\n", sep = "")
  cat("  did_status           : ", x$did_status, "\n", sep = "")
  if (!is.na(x$did_condition)) {
    cat("  did_condition        : ", x$did_condition, "\n", sep = "")
  }
  cat("  prior_canonical_kind : ", x$prior_canonical_kind, "\n", sep = "")
  invisible(x)
}

#' Construct a family object for AMM fitting
#'
#' Define the response distribution that links the individual parameter
#' theta_i to the observed outcome y_i. The family object carries the
#' link function, the inverse link, the metadata for the parameter
#' identifiability condition (D-ID) of Lemma 1B in Block 1, and the
#' family identifier consumed by the Stan code generator.
#'
#' @param name Character scalar identifying the family. One of
#'   \code{"gaussian"}, \code{"poisson"}, \code{"neg_binomial_2"},
#'   \code{"bernoulli"}, \code{"beta"}, \code{"gamma"},
#'   \code{"student_t"}, \code{"tweedie"}, \code{"zip"}, \code{"zinb"},
#'   \code{"hurdle_poisson"} or \code{"hurdle_neg_binomial_2"} for
#'   built-in families. For user-defined families, see
#'   \code{\link{gdpar_family_custom}}.
#' @param link Character scalar identifying the link function. The
#'   default is the canonical link for each family.
#' @param did_override Optional named list to override the canonical
#'   identifiability descriptors of one or more slots without changing
#'   the likelihood or links of the family. Keys are slot names of the
#'   family's \code{param_specs} (e.g., \code{"mu"}, \code{"sigma"},
#'   \code{"nu"} for Student-t); each value is a list with optional
#'   fields \code{did_status}, \code{did_condition},
#'   \code{did_reference}. \code{did_status} must be one of
#'   \code{"holds"}, \code{"holds_under_condition"},
#'   \code{"user_responsible"}. Defaults to \code{NULL} (canonical
#'   D-ID stays in place). Use this when your design violates the
#'   canonical D-ID assumption of the family and you want the pre-fit
#'   identifiability report to reflect that, without forking the
#'   likelihood through the custom-family surface.
#'
#' @return An object of class \code{gdpar_family} with components
#'   \code{name}, \code{link}, \code{inv_link}, \code{linkfun},
#'   \code{stan_id}, \code{has_dispersion}, \code{did_status},
#'   \code{did_condition} and \code{did_reference}. A \code{print}
#'   method provides a human-readable summary.
#'
#' @details
#' The Path 1 implementation in this version of the package fits the
#' AMM canonical form on the linear-predictor scale of the family. The
#' inverse link is applied at the likelihood block in Stan; the
#' centering, anchoring and prior specifications all live on the
#' linear-predictor scale, consistent with assumptions (C1)-(C6) of
#' Block 1 and the prior conditions (PRIOR-KL), (PRIOR-THICK) of Block 4.
#'
#' Built-in families and their D-ID status:
#' \describe{
#'   \item{\code{"gaussian"}}{Identity link by default. D-ID holds
#'     unconditionally for the location parameter when the variance is
#'     identifiable from the data.}
#'   \item{\code{"poisson"}}{Log link by default. D-ID holds
#'     unconditionally on the rate parameter.}
#'   \item{\code{"neg_binomial_2"}}{Log link by default; Stan's
#'     neg_binomial_2 parametrization with mean mu and dispersion phi
#'     such that variance equals mu + mu^2 / phi. D-ID holds
#'     unconditionally for mu when phi is identifiable from the data.}
#'   \item{\code{"bernoulli"}}{Logit link by default. D-ID holds
#'     unconditionally on the success probability.}
#'   \item{\code{"beta"}}{Logit link by default for the mean
#'     \code{mu} on (0, 1); precision \code{phi} on log link.
#'     Stan parametrization \code{beta_proportion(mu, phi)} with
#'     variance \code{mu*(1-mu)/(1+phi)}. D-ID holds for \code{mu}
#'     unconditionally; for \code{phi} when the empirical dispersion
#'     identifies it (sub-phase 8.3.4 of Block 8).}
#'   \item{\code{"gamma"}}{Log link by default for the mean
#'     \code{mu} on positive reals; shape on log link. Stan
#'     parametrization \code{gamma(shape, shape/mu)} with variance
#'     \code{mu^2/shape}. D-ID holds for \code{mu} unconditionally;
#'     for \code{shape} when the empirical dispersion identifies it
#'     (sub-phase 8.3.4 of Block 8).}
#'   \item{\code{"student_t"}}{Identity link by default for the
#'     location \code{mu} on the real line; \code{sigma} on log link;
#'     \code{nu} (degrees of freedom) on log link. Stan
#'     parametrization \code{student_t(nu, mu, sigma)} with mean
#'     \code{mu} (for \code{nu > 1}) and variance
#'     \code{sigma^2 * nu / (nu - 2)} (for \code{nu > 2}). D-ID holds
#'     for \code{mu} unconditionally; for \code{sigma} when the
#'     empirical residual dispersion identifies it; for \code{nu}
#'     when the empirical tail thickness identifies it. Sub-phase
#'     8.3.5a of Block 8 wires the family exclusively under the K = 3
#'     distributional regression path of \code{amm_distrib_K.stan};
#'     routings with K < 3 (population-scoped \code{sigma} and / or
#'     \code{nu}) are deferred.}
#'   \item{\code{"tweedie"}}{Log link by default for the mean
#'     \code{mu} on positive reals; dispersion \code{phi} on log link;
#'     power \code{p} on identity link with support on the open
#'     interval \code{(1.01, 1.99)} (canonical compound Poisson-gamma
#'     regime, with point mass at zero and continuous positive body).
#'     Stan parametrization \code{tweedie(mu, phi, p)} with mean
#'     \code{mu} and variance \code{phi * mu^p}. The log-pdf is not
#'     native in Stan and is implemented in the model's
#'     \code{functions} block as a hybrid: the Dunn--Smyth (2005)
#'     infinite-series expansion in the central region
#'     \code{|p - 1.5| < tau} (\code{tau = 0.4} by default) and the
#'     saddlepoint approximation elsewhere, with the same dispatch
#'     applied to the random-number generator. D-ID holds for
#'     \code{mu} unconditionally; for \code{phi} when the empirical
#'     dispersion identifies it; for \code{p} when the empirical
#'     balance between zero mass and continuous body identifies it.
#'     Sub-phase 8.3.5b of Block 8 wires the family exclusively under
#'     the K = 3 distributional regression path of
#'     \code{amm_distrib_K.stan}.}
#'   \item{\code{"zip"}}{Zero-inflated Poisson. Log link by default
#'     for the count mean \code{mu} on positive reals; mixture
#'     probability \code{pi} on logit link with support on
#'     \code{(0, 1)}. Likelihood: with probability \code{pi} the
#'     outcome is a structural zero, otherwise it is drawn from
#'     \code{Poisson(mu)}. Stan implementation via
#'     \code{log_sum_exp(bernoulli_logit_lpmf(1 | eta_pi),
#'     bernoulli_logit_lpmf(0 | eta_pi) + poisson_log_lpmf(0 | eta_mu))}
#'     for \code{y = 0} and \code{bernoulli_logit_lpmf(0 | eta_pi) +
#'     poisson_log_lpmf(y | eta_mu)} for \code{y > 0}. D-ID holds for
#'     \code{mu} unconditionally; for \code{pi} when the empirical
#'     proportion of structural zeros, distinct from the sampling
#'     zeros of the count component, is sufficient to identify the
#'     zero-inflation probability. Sub-phase 8.3.6 of Block 8 wires
#'     the family exclusively under the K = 2 distributional
#'     regression path of \code{amm_distrib_K.stan}; routings with
#'     K < 2 are deferred.}
#'   \item{\code{"zinb"}}{Zero-inflated negative binomial. Log link by
#'     default for the count mean \code{mu} on positive reals;
#'     dispersion \code{phi} on log link; mixture probability
#'     \code{pi} on logit link. Likelihood: with probability \code{pi}
#'     the outcome is a structural zero, otherwise it is drawn from
#'     \code{NegBinomial2(mu, phi)} with variance
#'     \code{mu + mu^2 / phi}. D-ID holds for \code{mu} unconditionally;
#'     for \code{phi} when the empirical overdispersion identifies it;
#'     for \code{pi} as in the zero-inflated Poisson case. Sub-phase
#'     8.3.6 wires the family exclusively under the K = 3
#'     distributional regression path; routings with K < 3 are
#'     deferred.}
#'   \item{\code{"hurdle_poisson"}}{Hurdle Poisson. Log link by
#'     default for the truncated-count mean \code{mu} on positive
#'     reals; hurdle probability \code{pi} on logit link. Likelihood:
#'     a Bernoulli draw with probability \code{pi} decides whether the
#'     outcome equals zero or is strictly positive; the positive
#'     branch is drawn from a Poisson truncated at one. Stan
#'     implementation via \code{bernoulli_logit_lpmf(1 | eta_pi)} for
#'     \code{y = 0} and \code{bernoulli_logit_lpmf(0 | eta_pi) +
#'     poisson_log_lpmf(y | eta_mu) - log1m_exp(-exp(eta_mu))} for
#'     \code{y > 0}. Distinct from the zero-inflated Poisson in that
#'     the zero mass is a structural decision, not a mixture between a
#'     structural zero and a Poisson sampling zero. D-ID holds for
#'     \code{mu} unconditionally; for \code{pi} when the empirical
#'     proportion of zeros is sufficient to identify the hurdle
#'     probability. Sub-phase 8.3.6 wires the family exclusively under
#'     the K = 2 distributional regression path.}
#'   \item{\code{"hurdle_neg_binomial_2"}}{Hurdle negative binomial.
#'     Log link by default for the truncated-count mean \code{mu};
#'     dispersion \code{phi} on log link; hurdle probability \code{pi}
#'     on logit link. Likelihood: Bernoulli draw decides zero vs.
#'     positive; the positive branch is drawn from a
#'     \code{NegBinomial2(mu, phi)} truncated at one. D-ID holds for
#'     \code{mu} unconditionally; for \code{phi} as in negative
#'     binomial; for \code{pi} as in hurdle Poisson. Sub-phase 8.3.6
#'     wires the family exclusively under the K = 3 distributional
#'     regression path.}
#' }
#'
#' @section Methodological notes:
#' Identifiability of the response family in its parameter is a
#' structural property of the model, not a property that can be tested
#' from finite data. The package therefore documents this property
#' through the \code{did_status} field but does not attempt to verify
#' it at fitting time. When the family carries a conditional D-ID
#' status, the fitting routine emits an informative message reminding
#' the user of the relevant condition. See Lemma 1B and the
#' commentary in Block 1, Section 6.4.
#'
#' The returned family declares every structural parameter the
#' distribution admits as \emph{eligible} for an individual
#' specification (e.g., \code{mu} and \code{sigma} for Gaussian;
#' \code{mu} and \code{phi} for negative binomial; \code{mu} for
#' Poisson and Bernoulli). Each eligible parameter carries a canonical
#' \code{scope} in its \code{gdpar_param_spec}: the location parameter
#' is \code{per_observation} (it is modeled through the AMM canonical
#' form) and auxiliary parameters default to \code{population}.
#' K-individual membership (which auxiliaries are promoted to
#' \code{per_observation}) is declared exclusively at the entry of
#' \code{\link{gdpar}} via a \code{\link{gdpar_formula_set}} (for the
#' high-level formula path) or via a named list of
#' \code{\link{amm_spec}} objects (for the low-level path). Parameters
#' that are not named in that entry keep their canonical scope and are
#' estimated as population-level constants. The family is the registry
#' of eligibles, and the entry to \code{gdpar()} is the single source
#' of truth for K.
#'
#' @section Dependencies:
#' This function uses \pkg{stats} family objects (\code{gaussian},
#' \code{poisson}, \code{binomial}) for the link metadata.
#'
#' @references
#' See \code{vignette("01_amm_identifiability", package = "gdpar")},
#' Section 6.4 (Lemma 1B) for D-ID; Lambert (1992), Greene (1994) and
#' Mullahy (1986) for identifiability of zero-inflated and hurdle
#' counts (sub-phase 8.3.6 of Block 8).
#'
#' @examples
#' fam <- gdpar_family("poisson")
#' print(fam)
#'
#' @seealso \code{\link{gdpar_family_custom}}, \code{\link{gdpar}}
#' @export
gdpar_family <- function(name = c("gaussian", "poisson",
                                  "neg_binomial_2", "bernoulli",
                                  "beta", "gamma", "student_t",
                                  "tweedie", "zip", "zinb",
                                  "hurdle_poisson",
                                  "hurdle_neg_binomial_2"),
                         link = NULL,
                         did_override = NULL) {
  name <- match.arg(name)
  default_link <- switch(name,
    gaussian              = "identity",
    poisson               = "log",
    neg_binomial_2        = "log",
    bernoulli             = "logit",
    beta                  = "logit",
    gamma                 = "log",
    student_t             = "identity",
    tweedie               = "log",
    zip                   = "log",
    zinb                  = "log",
    hurdle_poisson        = "log",
    hurdle_neg_binomial_2 = "log"
  )
  if (is.null(link)) {
    link <- default_link
  }
  allowed_links <- list(
    gaussian              = c("identity", "log"),
    poisson               = "log",
    neg_binomial_2        = "log",
    bernoulli             = "logit",
    beta                  = "logit",
    gamma                 = "log",
    student_t             = c("identity", "log"),
    tweedie               = "log",
    zip                   = "log",
    zinb                  = "log",
    hurdle_poisson        = "log",
    hurdle_neg_binomial_2 = "log"
  )
  if (!link %in% allowed_links[[name]]) {
    gdpar_abort(
      sprintf(
        "Family '%s' does not support link '%s'. Allowed links: %s.",
        name, link, paste(sQuote(allowed_links[[name]]), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(family = name, link = link,
                  allowed = allowed_links[[name]])
    )
  }
  param_specs <- .gdpar_family_param_specs_for(name, link)
  param_specs <- .gdpar_apply_did_override(param_specs, did_override, name)
  stan_id <- switch(name,
    gaussian              = 1L,
    poisson               = 2L,
    neg_binomial_2        = 3L,
    bernoulli             = 4L,
    beta                  = 5L,
    gamma                 = 6L,
    student_t             = 8L,
    tweedie               = 9L,
    zip                   = 10L,
    zinb                  = 11L,
    hurdle_poisson        = 12L,
    hurdle_neg_binomial_2 = 13L
  )
  has_dispersion <- length(param_specs) > 1L
  first <- param_specs[[1L]]
  obj <- list(
    name           = name,
    link           = first$link,
    inv_link       = first$inv_link,
    linkfun        = first$linkfun,
    stan_id        = stan_id,
    has_dispersion = has_dispersion,
    did_status     = first$did_status,
    did_condition  = first$did_condition,
    did_reference  = first$did_reference,
    param_specs    = param_specs
  )
  class(obj) <- c("gdpar_family", "list")
  obj
}

#' Construct a custom family object for AMM fitting
#'
#' Build a user-defined family for use with \code{\link{gdpar}} when
#' the built-in families of \code{\link{gdpar_family}} do not cover the
#' application. The user is responsible for declaring whether the
#' identifiability condition (D-ID) of Lemma 1B in Block 1 holds for
#' the family.
#'
#' @param name Character scalar identifying the custom family. Must
#'   not coincide with any built-in family name.
#' @param link Character scalar identifying the link function. One of
#'   \code{"identity"}, \code{"log"}, \code{"logit"}.
#' @param did_holds Logical scalar. The user must explicitly declare
#'   whether the family is identifiable in its parameter; a missing
#'   declaration raises an error.
#' @param did_condition Character scalar describing any condition under
#'   which D-ID holds when \code{did_holds = TRUE} but identifiability
#'   is conditional. Use \code{NA_character_} if D-ID holds
#'   unconditionally.
#' @param stan_loglik_block Character scalar with a Stan code snippet
#'   for the \code{model} block: declares
#'   \code{target += custom_lpdf | custom_lpmf} for one observation.
#'   The snippet must reference the linear predictor \code{eta[i]} and
#'   either \code{y_real[i]} (when \code{y_type = "real"}) or
#'   \code{y_int[i]} (when \code{y_type = "integer"}); the legacy
#'   placeholder \code{y[i]} is rejected.
#' @param stan_log_lik_block Character scalar with the
#'   \code{generated quantities} snippet that assigns to
#'   \code{log_lik[i]} for one observation, e.g.
#'   \code{log_lik[i] = normal_lpdf(log(y_real[i]) | eta[i], sigma);}.
#'   Used by \code{\link{gdpar_loo}} downstream.
#' @param stan_y_pred_block Character scalar with the
#'   \code{generated quantities} snippet that assigns to
#'   \code{y_pred[i]} for one observation, e.g.
#'   \code{y_pred[i] = exp(normal_rng(eta[i], sigma));}. Used by
#'   posterior-predictive utilities.
#' @param y_type Character scalar, one of \code{"real"} or
#'   \code{"integer"}, declaring whether the outcome is real-valued
#'   (Stan template references \code{y_real[i]}) or integer-valued
#'   (Stan template references \code{y_int[i]}).
#' @param did_reference Character scalar with a citation supporting the
#'   D-ID declaration.
#'
#' @return An object of class \code{gdpar_family}.
#'
#' @details
#' Building a custom family is an advanced use of the package. The user
#' assumes the responsibility of ensuring (i) that the Stan likelihood
#' block is mathematically correct and (ii) that the family is
#' identifiable in its parameter. The package emits an informative
#' message when the custom family is created, restating these
#' responsibilities.
#'
#' @section Methodological notes:
#' The package never attempts to test identifiability from data; it
#' only registers the user's declaration. See Lemma 1B in Block 1,
#' Section 6.4.
#'
#' @section Dependencies:
#' None beyond the \pkg{base} R installation.
#'
#' @examples
#' my_family <- gdpar_family_custom(
#'   name               = "my_log_normal",
#'   link               = "log",
#'   did_holds          = TRUE,
#'   did_condition      = NA_character_,
#'   stan_loglik_block  =
#'     "target += normal_lpdf(log(y_real[i]) | eta[i], sigma_y[1]);",
#'   stan_log_lik_block =
#'     "log_lik[i] = normal_lpdf(log(y_real[i]) | eta[i], sigma_y[1]);",
#'   stan_y_pred_block  =
#'     "y_pred[i] = exp(normal_rng(eta[i], sigma_y[1]));",
#'   y_type             = "real",
#'   did_reference      = "User declaration"
#' )
#' print(my_family)
#'
#' @seealso \code{\link{gdpar_family}}
#' @export
gdpar_family_custom <- function(name, link, did_holds, did_condition,
                                stan_loglik_block,
                                stan_log_lik_block,
                                stan_y_pred_block,
                                y_type,
                                did_reference) {
  if (!is.character(name) || length(name) != 1L) {
    gdpar_abort("Argument 'name' must be a character scalar.",
                class = "gdpar_input_error")
  }
  builtins <- c("gaussian", "poisson", "neg_binomial_2", "bernoulli",
                "beta", "gamma", "student_t", "tweedie",
                "zip", "zinb", "hurdle_poisson",
                "hurdle_neg_binomial_2")
  if (name %in% builtins) {
    gdpar_abort(
      sprintf(
        "Custom family name must not coincide with a built-in family ('%s').",
        name
      ),
      class = "gdpar_input_error"
    )
  }
  if (!is.character(link) || !link %in% c("identity", "log", "logit")) {
    gdpar_abort(
      "Argument 'link' must be one of 'identity', 'log', 'logit'.",
      class = "gdpar_input_error"
    )
  }
  if (!is.logical(did_holds) || length(did_holds) != 1L ||
      is.na(did_holds)) {
    gdpar_abort(
      "Argument 'did_holds' must be TRUE or FALSE; missing is not allowed.",
      class = "gdpar_input_error"
    )
  }
  if (!is.character(did_condition) || length(did_condition) != 1L) {
    gdpar_abort(
      "Argument 'did_condition' must be a character scalar (NA_character_ if unconditional).",
      class = "gdpar_input_error"
    )
  }
  for (nm in c("stan_loglik_block", "stan_log_lik_block",
               "stan_y_pred_block")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1L || nchar(val) == 0L) {
      gdpar_abort(
        sprintf(
          "Argument '%s' must be a non-empty character scalar (Stan code snippet).",
          nm
        ),
        class = "gdpar_input_error",
        data = list(argument = nm)
      )
    }
  }
  if (!is.character(y_type) || length(y_type) != 1L ||
      !y_type %in% c("real", "integer")) {
    gdpar_abort(
      "Argument 'y_type' must be either \"real\" or \"integer\".",
      class = "gdpar_input_error",
      data = list(received = y_type)
    )
  }
  legacy_y_pattern <- "(^|[^_[:alnum:]])y\\["
  for (nm in c("stan_loglik_block", "stan_log_lik_block",
               "stan_y_pred_block")) {
    val <- get(nm)
    if (grepl(legacy_y_pattern, val)) {
      gdpar_abort(
        sprintf(
          paste(
            "Argument '%s' uses the legacy placeholder 'y[i]'. The Stan",
            "template references the outcome as 'y_real[i]' when",
            "y_type = \"real\" or 'y_int[i]' when y_type = \"integer\";",
            "update the snippet accordingly."
          ),
          nm
        ),
        class = "gdpar_input_error",
        data = list(argument = nm)
      )
    }
  }
  did_status <- if (did_holds && is.na(did_condition)) "holds"
                else if (did_holds) "holds_under_condition"
                else "user_responsible"
  custom_spec <- .gdpar_param_spec(
    name                 = "mu",
    link                 = link,
    did_status           = did_status,
    did_condition        = did_condition,
    did_reference        = did_reference,
    support              = "custom",
    prior_canonical_kind = "mu",
    scope                = "per_observation",
    family_role          = "location"
  )
  obj <- list(
    name               = name,
    link               = custom_spec$link,
    inv_link           = custom_spec$inv_link,
    linkfun            = custom_spec$linkfun,
    stan_id            = NA_integer_,
    has_dispersion     = FALSE,
    did_status         = custom_spec$did_status,
    did_condition      = custom_spec$did_condition,
    did_reference      = custom_spec$did_reference,
    stan_loglik_block  = stan_loglik_block,
    stan_log_lik_block = stan_log_lik_block,
    stan_y_pred_block  = stan_y_pred_block,
    y_type             = y_type,
    is_custom          = TRUE,
    param_specs        = list(custom_spec)
  )
  class(obj) <- c("gdpar_family", "list")
  gdpar_inform(
    sprintf(
      paste(
        "Custom family '%s' registered. The user is responsible for",
        "the correctness of the Stan likelihood block and for the",
        "identifiability declaration. The package documents this",
        "declaration but does not verify it from data."
      ),
      name
    ),
    class = "gdpar_did_message"
  )
  obj
}

#' Internal: registry of canonical K = 2 custom-family patterns
#'
#' Returns a named list of canonical bi-parametric (K = 2) likelihood
#' patterns that are wired into \code{inst/stan/amm_distrib_K.stan} and
#' that can be exposed under a user-chosen family name through
#' \code{\link{gdpar_family_custom_K}}. Each entry carries the
#' \code{stan_id} branch already present in the Stan template, the
#' canonical \code{param_specs} of the two slots (with their
#' \code{link} and \code{prior_canonical_kind}), the outcome storage
#' type (\code{y_type}), and the canonical D-ID status.
#'
#' Sub-phase 8.3.4 of Block 8 opens this registry with a single new
#' pattern (\code{lognormal_loc_scale}) beyond the four built-in
#' distributional regression branches. The registry is the canonical
#' surface through which the descriptor-based custom-family wiring
#' (D-A3.B; option (b) of the sub-phase scoping) operates: a user
#' selects a \code{stan_lpdf_id} from this whitelist and the
#' constructor instantiates a \code{gdpar_family} with the registered
#' slot configuration and the \code{stan_id} of the corresponding
#' template branch. No free-form Stan code is accepted.
#'
#' Future sub-phases extend the registry by (i) adding a new
#' \code{stan_id} branch to \code{amm_distrib_K.stan}, (ii) updating
#' the whitelist in \code{.assemble_stan_data_K()} and
#' \code{.gdpar_guard_K_below_family_min()}, and (iii) appending a
#' row here.
#'
#' @return Named list (key = \code{stan_lpdf_id}) of pattern
#'   descriptors with components \code{stan_id} (integer),
#'   \code{slot_specs} (length-2 list of named lists with
#'   \code{name}, \code{link}, \code{support},
#'   \code{prior_canonical_kind}, \code{scope}, \code{family_role},
#'   \code{did_status}, \code{did_condition},
#'   \code{did_reference}), \code{y_type} (character).
#' @keywords internal
#' @noRd
.gdpar_K_custom_patterns <- function() {
  list(
    lognormal_loc_scale = list(
      stan_id = 7L,
      y_type  = "real",
      slot_specs = list(
        list(
          name                 = "mu",
          link                 = "identity",
          support              = "real_line",
          prior_canonical_kind = "mu",
          scope                = "per_observation",
          family_role          = "location",
          did_status           = "holds",
          did_condition        = NA_character_,
          did_reference        = "Block 1, Section 6.4 (Lemma 1B)"
        ),
        list(
          name                 = "sigma",
          link                 = "log",
          support              = "positive_real",
          prior_canonical_kind = "log_sigma",
          scope                = "population",
          family_role          = "scale",
          did_status           = "holds_under_condition",
          did_condition        = "empirical log-scale dispersion sufficient to identify sigma",
          did_reference        = "Block 1, Section 6.4 (Lemma 1B)"
        )
      )
    )
  )
}

#' Construct a K = 2 custom family from a canonical lpdf pattern
#'
#' Build a custom distributional regression family by selecting a
#' canonical bi-parametric likelihood pattern from the registry
#' returned by \code{.gdpar_K_custom_patterns}. The constructor
#' is the descriptor-based wiring agreed in sub-phase 8.3.4 of Block 8
#' (D-A3.B; option (b) of the scoping): the user does \emph{not}
#' supply free-form Stan code; they choose a \code{stan_lpdf_id} from
#' the whitelist and the family routes through the same
#' \code{amm_distrib_K.stan} branch that the built-in distributional
#' regression families use.
#'
#' This is the K = 2 sibling of \code{\link{gdpar_family_custom}} (the
#' K = 1 free-form custom path). The two coexist with distinct
#' contracts: \code{gdpar_family_custom()} accepts user-authored Stan
#' code blocks at K = 1 (user-validated identifiability and
#' correctness); \code{gdpar_family_custom_K()} accepts only registered
#' canonical patterns at K = 2 (descriptor-validated, no Stan
#' injection surface).
#'
#' @param name Character scalar identifying the custom family. Must
#'   not coincide with any built-in family name or another registered
#'   custom-K family in the calling session.
#' @param stan_lpdf_id Character scalar selecting one of the canonical
#'   patterns in \code{.gdpar_K_custom_patterns}. Sub-phase
#'   8.3.4 opens the registry with \code{"lognormal_loc_scale"}.
#' @param did_holds Logical scalar declaring whether the D-ID
#'   condition of Lemma 1B holds for the family. Default \code{TRUE}
#'   (the registry's pattern-level declaration is used; see Details).
#' @param did_condition Optional character scalar to override the
#'   pattern-level D-ID condition. Default \code{NULL} (use registry).
#' @param did_reference Optional character scalar with the user's
#'   citation supporting the D-ID declaration. Default \code{NULL}
#'   (use registry).
#'
#' @return An object of class \code{gdpar_family} with two
#'   \code{param_specs} (slot 1 = location, slot 2 = scale), the
#'   \code{stan_id} of the pattern (e.g. 7 for
#'   \code{lognormal_loc_scale}), and \code{is_custom = TRUE}. The
#'   family routes through \code{amm_distrib_K.stan} via the
#'   \code{family_id_k} dispatcher.
#'
#' @details
#' The descriptor path enforces structural validation: \code{name}
#' must be unique, \code{stan_lpdf_id} must be in the registry, and
#' the slot configuration (links, supports, prior kinds) is fixed by
#' the registry entry. Users who need to deviate from the canonical
#' parametrization must request a new entry in the registry; this
#' keeps the Stan-side surface auditable and bit-exact across
#' versions of the package.
#'
#' Sub-phase 8.3.4 (D-A3.B option (b)) decided this descriptor
#' approach over (i) plug-in of free-form Stan code per family slot
#' and (ii) a hybrid descriptor + escape-hatch. The decision applies
#' \code{[[feedback-max-robustness-priority]]} (structural validation
#' over surface flexibility) at the cost of forcing the user to
#' contribute upstream when a new K = 2 family is needed.
#'
#' @section Methodological notes:
#' The K = 2 custom family inherits the canonical likelihood, links,
#' and priors of the chosen registry pattern. The user only chooses
#' the family name and (optionally) overrides the D-ID metadata. The
#' AMM canonical form, the per-slot \code{scope}, and the dispatch in
#' \code{amm_distrib_K.stan} are the same as for the built-in
#' \code{gdpar_family} families wired in 8.3.4.
#'
#' @examples
#' my_lognorm <- gdpar_family_custom_K(
#'   name          = "my_lognormal_K2",
#'   stan_lpdf_id  = "lognormal_loc_scale",
#'   did_holds     = TRUE,
#'   did_reference = "User declaration"
#' )
#' print(my_lognorm)
#'
#' @seealso \code{\link{gdpar_family}}, \code{\link{gdpar_family_custom}}
#' @export
gdpar_family_custom_K <- function(name, stan_lpdf_id,
                                   did_holds = TRUE,
                                   did_condition = NULL,
                                   did_reference = NULL) {
  if (!is.character(name) || length(name) != 1L || nchar(name) == 0L) {
    gdpar_abort("Argument 'name' must be a non-empty character scalar.",
                class = "gdpar_input_error")
  }
  builtins <- c("gaussian", "poisson", "neg_binomial_2", "bernoulli",
                "beta", "gamma", "student_t", "tweedie",
                "zip", "zinb", "hurdle_poisson",
                "hurdle_neg_binomial_2")
  if (name %in% builtins) {
    gdpar_abort(
      sprintf(
        "Custom family name must not coincide with a built-in family ('%s').",
        name
      ),
      class = "gdpar_input_error",
      data = list(name = name)
    )
  }
  if (!is.character(stan_lpdf_id) || length(stan_lpdf_id) != 1L) {
    gdpar_abort(
      "Argument 'stan_lpdf_id' must be a character scalar.",
      class = "gdpar_input_error"
    )
  }
  patterns <- .gdpar_K_custom_patterns()
  if (!stan_lpdf_id %in% names(patterns)) {
    gdpar_abort(
      sprintf(
        paste(
          "stan_lpdf_id '%s' is not registered. Allowed patterns: %s.",
          "Sub-phase 8.3.4 of Block 8 opens the registry with",
          "'lognormal_loc_scale'; future sub-phases extend the whitelist."
        ),
        stan_lpdf_id,
        paste(sQuote(names(patterns)), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(stan_lpdf_id = stan_lpdf_id,
                  allowed = names(patterns))
    )
  }
  if (!is.logical(did_holds) || length(did_holds) != 1L ||
      is.na(did_holds)) {
    gdpar_abort(
      "Argument 'did_holds' must be TRUE or FALSE; missing is not allowed.",
      class = "gdpar_input_error"
    )
  }
  pattern <- patterns[[stan_lpdf_id]]
  param_specs <- lapply(pattern$slot_specs, function(s) {
    did_status_s <- if (!did_holds) "user_responsible" else s$did_status
    did_cond_s <- if (!is.null(did_condition)) did_condition else s$did_condition
    did_ref_s <- if (!is.null(did_reference)) did_reference else s$did_reference
    .gdpar_param_spec(
      name                 = s$name,
      link                 = s$link,
      did_status           = did_status_s,
      did_condition        = did_cond_s,
      did_reference        = did_ref_s,
      support              = s$support,
      prior_canonical_kind = s$prior_canonical_kind,
      scope                = s$scope,
      family_role          = s$family_role
    )
  })
  first <- param_specs[[1L]]
  obj <- list(
    name           = name,
    link           = first$link,
    inv_link       = first$inv_link,
    linkfun        = first$linkfun,
    stan_id        = as.integer(pattern$stan_id),
    has_dispersion = TRUE,
    did_status     = first$did_status,
    did_condition  = first$did_condition,
    did_reference  = first$did_reference,
    y_type         = pattern$y_type,
    is_custom      = TRUE,
    stan_lpdf_id   = stan_lpdf_id,
    param_specs    = param_specs
  )
  class(obj) <- c("gdpar_family", "list")
  gdpar_inform(
    sprintf(
      paste(
        "Custom K = 2 family '%s' registered using canonical pattern '%s'.",
        "The Stan-side wiring is bit-identical to the registry entry;",
        "identifiability follows the pattern declaration unless",
        "explicitly overridden by the user."
      ),
      name, stan_lpdf_id
    ),
    class = "gdpar_did_message"
  )
  obj
}

#' Print method for gdpar_family objects
#'
#' @param x An object of class \code{gdpar_family}.
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_family <- function(x, ...) {
  cat("<gdpar_family>\n")
  cat("  name           : ", x$name, "\n", sep = "")
  cat("  link           : ", x$link, "\n", sep = "")
  cat("  has_dispersion : ", x$has_dispersion, "\n", sep = "")
  cat("  did_status     : ", x$did_status, "\n", sep = "")
  if (!is.na(x$did_condition)) {
    cat("  did_condition  : ", x$did_condition, "\n", sep = "")
  }
  cat("  did_reference  : ", x$did_reference, "\n", sep = "")
  if (!is.null(x$param_specs)) {
    spec_names <- vapply(x$param_specs, function(s) s$name, character(1L))
    spec_scopes <- vapply(x$param_specs, function(s) s$scope, character(1L))
    spec_summary <- paste(
      mapply(function(n, s) sprintf("%s (%s)", n, s),
             spec_names, spec_scopes),
      collapse = ", "
    )
    cat("  param_specs    : ", spec_summary, "\n", sep = "")
  }
  invisible(x)
}

#' Construct a multivariate family for AMM fitting with p > 1
#'
#' Build a per-coordinate family object for use with \code{\link{gdpar}}
#' when the AMM specification has dimension \code{p > 1L}. The resulting
#' object declares one univariate family per coordinate of the
#' individual parameter vector \eqn{\theta_i \in \mathbb{R}^p}; the
#' likelihood factorizes across coordinates as
#' \deqn{p(y_i \mid \theta_i) = \prod_{k=1}^{p} D_k(y_{ik} \mid \theta_i[k]),}
#' with cross-dimensional coupling carried exclusively by the modulating
#' component \eqn{W(\theta_{\mathrm{ref}})} of the AMM canonical form.
#'
#' @param family Either a character scalar with the name of a built-in
#'   family (one of \code{"gaussian"}, \code{"poisson"},
#'   \code{"neg_binomial_2"}, \code{"bernoulli"}), an object of class
#'   \code{gdpar_family} produced by \code{\link{gdpar_family}}, or a
#'   list of length \code{p} where each entry is itself a
#'   \code{gdpar_family} (heterogeneous coordinates). In this version
#'   the list form is restricted to the homogeneous case (all entries
#'   must share the same \code{stan_id}); heterogeneous families per
#'   coordinate are deferred to a later sub-phase.
#' @param p Positive integer giving the dimension of \eqn{\theta_i}.
#'   Must match the \code{p} of the \code{\link{amm_spec}} passed to
#'   \code{\link{gdpar}}.
#' @param link Character scalar identifying the link function when
#'   \code{family} is supplied as a name. Ignored when \code{family}
#'   is already a \code{gdpar_family} object or a list. Defaults to the
#'   canonical link of the named family.
#'
#' @return An object of class \code{gdpar_family_multi} with components
#'   \code{families} (a list of \code{p} \code{gdpar_family} objects),
#'   \code{p}, \code{homogeneous}, \code{stan_id} (the common Stan
#'   family identifier when homogeneous), \code{has_dispersion} (the
#'   common dispersion flag), \code{name} (the common family name) and
#'   \code{did_status} (the common identifiability status). A
#'   \code{print} method provides a human-readable summary.
#'
#' @details
#' The factorization is the canonization of architectural Option B of
#' the Phase F decision (handoff 10 continuation, 2026-05-11):
#' \eqn{y_i \in \mathbb{R}^p} is multivariate with marginals independent
#' conditional on \eqn{\theta_i}; cross-dimensional dependence enters
#' the model through the coupling of \eqn{\theta_i[k]} via
#' \eqn{W(\theta_{\mathrm{ref}})}. Multi-parametric families (a single
#' univariate outcome parametrized by the whole vector
#' \eqn{\theta_i \in \mathbb{R}^p}, e.g., gaussian with
#' \eqn{\theta_i = (\mu_i, \log \sigma_i)} in the distributional
#' regression sense) are deferred to a dedicated post-validation block;
#' see the project memory entry
#' \code{project_gdpar_multiparametric_extension_postvalidation}.
#'
#' @section Methodological notes:
#' The identifiability condition (D-ID) of Lemma 1B in Block 1 applies
#' coordinate-wise under this factorization: each univariate marginal
#' \eqn{D_k} identifies \eqn{\theta_i[k]} from \eqn{y_{ik}}
#' independently. The cross-dimensional identifiability condition
#' (C4-bis), which guards against aliasing between coordinates of
#' \eqn{\theta_{\mathrm{ref}}} that share basis structure, is checked
#' by \code{\link{gdpar_check_identifiability}} (Phase H pending).
#'
#' @section Dependencies:
#' Calls \code{\link{gdpar_family}} when \code{family} is supplied as
#' a name.
#'
#' @examples
#' fam_mv <- gdpar_family_multi("gaussian", p = 2L)
#' print(fam_mv)
#'
#' fam_mv2 <- gdpar_family_multi(gdpar_family("poisson"), p = 3L)
#' print(fam_mv2)
#'
#' @seealso \code{\link{gdpar_family}}, \code{\link{amm_spec}}
#' @export
gdpar_family_multi <- function(family, p, link = NULL) {
  assert_count(p, "p")
  p <- as.integer(p)
  if (p < 1L) {
    gdpar_abort(
      "Argument 'p' must be a positive integer.",
      class = "gdpar_input_error",
      data = list(p = p)
    )
  }
  if (is.character(family) && length(family) == 1L) {
    base_family <- gdpar_family(name = family, link = link)
    families <- replicate(p, base_family, simplify = FALSE)
  } else if (inherits(family, "gdpar_family")) {
    if (!is.null(link)) {
      gdpar_abort(
        paste(
          "Argument 'link' is ignored when 'family' is a gdpar_family",
          "object; the link of the supplied object is used. Supply",
          "'family' as a name string to specify a non-default link."
        ),
        class = "gdpar_input_error"
      )
    }
    families <- replicate(p, family, simplify = FALSE)
  } else if (is.list(family)) {
    if (length(family) != p) {
      gdpar_abort(
        sprintf(
          paste(
            "When 'family' is a list, its length must equal p = %d;",
            "received length %d."
          ),
          p, length(family)
        ),
        class = "gdpar_input_error",
        data = list(p = p, received_length = length(family))
      )
    }
    for (k in seq_len(p)) {
      if (!inherits(family[[k]], "gdpar_family")) {
        gdpar_abort(
          sprintf(
            paste(
              "Entry %d of 'family' list must be a gdpar_family object;",
              "received an object of class %s."
            ),
            k, paste(sQuote(class(family[[k]])), collapse = ", ")
          ),
          class = "gdpar_input_error"
        )
      }
    }
    families <- family
  } else {
    gdpar_abort(
      paste(
        "Argument 'family' must be a name string (built-in family),",
        "a gdpar_family object, or a list of p gdpar_family objects."
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(family))
    )
  }
  stan_ids <- vapply(families, function(f) f$stan_id, integer(1L))
  homogeneous <- length(unique(stan_ids)) == 1L
  if (!homogeneous) {
    gdpar_abort(
      paste(
        "Heterogeneous families per coordinate are deferred to a later",
        "sub-phase of Block 5.2 (Phase F.x). All entries must currently",
        "share the same Stan family identifier. Received stan_ids: ",
        paste(stan_ids, collapse = ", "), "."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(stan_ids = stan_ids)
    )
  }
  links <- vapply(families, function(f) f$link, character(1L))
  if (length(unique(links)) > 1L) {
    gdpar_abort(
      paste(
        "All coordinates must currently share the same link function.",
        "Received links: ", paste(unique(links), collapse = ", "), "."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(links = links)
    )
  }
  param_specs_per_coord <- lapply(families, function(f) f$param_specs)
  obj <- list(
    families              = families,
    p                     = p,
    homogeneous           = homogeneous,
    stan_id               = families[[1L]]$stan_id,
    name                  = families[[1L]]$name,
    link                  = families[[1L]]$link,
    has_dispersion        = families[[1L]]$has_dispersion,
    did_status            = families[[1L]]$did_status,
    did_condition         = families[[1L]]$did_condition,
    did_reference         = families[[1L]]$did_reference,
    param_specs_per_coord = param_specs_per_coord
  )
  class(obj) <- c("gdpar_family_multi", "list")
  obj
}

#' Print method for gdpar_family_multi objects
#'
#' @param x An object of class \code{gdpar_family_multi}.
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_family_multi <- function(x, ...) {
  cat("<gdpar_family_multi> coord-wise factorization\n")
  cat("  p              : ", x$p, "\n", sep = "")
  cat("  homogeneous    : ", x$homogeneous, "\n", sep = "")
  cat("  name           : ", x$name, "\n", sep = "")
  cat("  link           : ", x$link, "\n", sep = "")
  cat("  has_dispersion : ", x$has_dispersion, "\n", sep = "")
  cat("  did_status     : ", x$did_status, "\n", sep = "")
  if (!is.na(x$did_condition)) {
    cat("  did_condition  : ", x$did_condition, "\n", sep = "")
  }
  cat("  did_reference  : ", x$did_reference, "\n", sep = "")
  if (!is.null(x$param_specs_per_coord) && x$homogeneous) {
    first_specs <- x$param_specs_per_coord[[1L]]
    spec_names <- vapply(first_specs, function(s) s$name, character(1L))
    spec_scopes <- vapply(first_specs, function(s) s$scope, character(1L))
    spec_summary <- paste(
      mapply(function(n, s) sprintf("%s (%s)", n, s),
             spec_names, spec_scopes),
      collapse = ", "
    )
    cat("  param_specs    : ", spec_summary,
        " per coord (x", x$p, ")\n", sep = "")
  }
  invisible(x)
}

#' Internal: promote a set of family parameters to per-observation scope
#'
#' Returns a copy of the supplied family in which every \code{param_spec}
#' whose \code{name} appears in \code{k_names} has its \code{scope}
#' field set to \code{"per_observation"}. Parameters whose name is not
#' in \code{k_names} retain their canonical \code{scope}. The function
#' is the K-individual promotion mechanism canonized in decision D of
#' sub-phase 8.3.3 (memory entry
#' \code{project_gdpar_block_8_3_extended_plan}, 2026-05-20):
#' \code{gdpar_family()} declares the set of eligible parameters; the
#' entry to \code{\link{gdpar}} (a \code{\link{gdpar_formula_set}} or a
#' named list of \code{\link{amm_spec}}) declares which subset is
#' individual, and this helper materializes the scope promotion on a
#' family copy that \code{\link{gdpar}} forwards to the codegen. The
#' helper \code{.gdpar_n_params_individual()} consumes the promoted
#' scope to compute K downstream.
#'
#' For \code{gdpar_family_multi} the promotion is applied to every
#' per-coordinate family in the \code{families} list (homogeneous
#' coordinates share the same eligible parameter set in this version).
#'
#' Validation: \code{k_names} must be a subset of the names of the
#' family's \code{param_specs}. A mismatch aborts with
#' \code{gdpar_input_error} naming the offending parameter and listing
#' the eligible names. An empty or \code{NULL} \code{k_names} returns
#' the family unchanged.
#'
#' @param family A \code{gdpar_family} or \code{gdpar_family_multi}.
#' @param k_names Character vector with the names of the parameters
#'   to promote. May be empty or \code{NULL} (no-op).
#'
#' @return A family object of the same class as \code{family} with the
#'   promoted \code{scope} fields. The original input is not mutated.
#'
#' @keywords internal
#' @noRd
.gdpar_promote_scope_per_observation <- function(family, k_names) {
  if (is.null(k_names) || length(k_names) == 0L) {
    return(family)
  }
  if (!is.character(k_names)) {
    gdpar_abort(
      "Argument 'k_names' must be a character vector.",
      class = "gdpar_internal_error",
      data = list(received = class(k_names))
    )
  }
  promote_one <- function(fam) {
    if (is.null(fam$param_specs)) {
      gdpar_abort(
        paste(
          "Internal error: family object has no 'param_specs' field;",
          "cannot promote scope."
        ),
        class = "gdpar_internal_error"
      )
    }
    eligible <- vapply(fam$param_specs, function(s) s$name, character(1L))
    bad <- setdiff(k_names, eligible)
    if (length(bad) > 0L) {
      gdpar_abort(
        sprintf(
          paste(
            "Parameter(s) %s named for K-individual specification are not",
            "eligible in family '%s'. Eligible parameters: %s."
          ),
          paste(sQuote(bad), collapse = ", "),
          fam$name %||% "<unknown>",
          paste(sQuote(eligible), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(unknown = bad, eligible = eligible,
                    family = fam$name)
      )
    }
    for (i in seq_along(fam$param_specs)) {
      if (fam$param_specs[[i]]$name %in% k_names) {
        fam$param_specs[[i]]$scope <- "per_observation"
      }
    }
    fam
  }
  if (inherits(family, "gdpar_family_multi")) {
    family$families <- lapply(family$families, promote_one)
    if (!is.null(family$param_specs_per_coord)) {
      family$param_specs_per_coord <- lapply(
        family$families, function(f) f$param_specs
      )
    }
    return(family)
  }
  promote_one(family)
}

# -----------------------------------------------------------------
# Sub-phase 8.3.7 (2026-05-21): heterogeneous families per slot.
#
# Decisions D1-D7 of the scoping (see project memory entry
# project_gdpar_block_8_3_7_closed; this comment block summarizes the
# operational shape):
#
# - D1 (A): family per slot only; the link of each slot stays canonical
#   for the family that governs the slot.
# - D2 (alpha): the family of the slot determines the parametrization
#   of the slot (link + canonical prior + support + did_status); the
#   likelihood is governed exclusively by the family of the location
#   slot. There is one observed y per fit.
# - D3 (I): Stan dispatch happens slot-by-slot via the data field
#   inv_link_id_per_slot, populated R-side.
# - D3.5 (L1 refined): when family_id_k[k] != family_id_k[1] (genuinely
#   heterogeneous), the slot inherits the link + canonical prior +
#   support + did_status of slot 1 (location) of the heterogeneous
#   family. When family_id_k[k] == family_id_k[1] (homogeneous), the
#   slot keeps the canonical link of slot k of the location family.
#   This refinement preserves the homogeneous regime bit-for-bit while
#   enabling cross-slot heterogeneity in the genuinely heterogeneous
#   case.
# - D4 (a)+(ii): strict validator on param_spec$support;
#   gdpar_input_error with rich messages that cite slot role, emitted
#   support, required support, and alternative families that satisfy
#   the predicate.
# - D5 (a): public API accepts a named list in the family argument;
#   names must match the slot names of the gdpar_formula_set.
# - D6 (a2): three smokes gated by env vars (Gauss+Beta scale /
#   Gauss+Gamma scale / NB+Beta dispersion).
# - D7 (a): clean refactor of the likelihood branches via the helper
#   apply_inv_link_by_id; goldens K=2 Beta/Gamma re-bootstrap at close.

#' Internal: canonical inverse-link ID per built-in family slot 1
#'
#' Returns an integer ID identifying the canonical inverse link of the
#' location slot (slot 1) of the family with the supplied stan_id.
#' The ID is consumed by the Stan helper apply_inv_link_by_id() inside
#' amm_distrib_K.stan: 0 = identity, 1 = inv_logit, 2 = exp. See D-D3
#' (decision 3) and D-D3.5 (L1 refined) of the 8.3.7 scoping.
#'
#' @param stan_id Integer scalar; one of 1, 3, 5, 6, 7, 8, 9, 10, 11,
#'   12, 13 corresponding to built-in K-individual families. Custom
#'   families with stan_id < 0 or NA fall back to a NULL error.
#' @return Integer scalar in {0, 1, 2}.
#' @keywords internal
#' @noRd
.gdpar_canonical_inv_link_id_slot1 <- function(stan_id) {
  if (is.null(stan_id) || is.na(stan_id)) {
    gdpar_abort(
      "Internal error: canonical inv_link id requested for NA stan_id.",
      class = "gdpar_internal_error"
    )
  }
  switch(as.character(as.integer(stan_id)),
    "1"  = 0L,   # gaussian:              mu identity
    "3"  = 2L,   # neg_binomial_2:        mu log
    "5"  = 1L,   # beta:                  mu logit
    "6"  = 2L,   # gamma:                 mu log
    "7"  = 0L,   # lognormal_loc_scale:   mu identity (on log y)
    "8"  = 0L,   # student_t:             mu identity
    "9"  = 2L,   # tweedie:               mu log
    "10" = 2L,   # zip:                   mu log
    "11" = 2L,   # zinb:                  mu log
    "12" = 2L,   # hurdle_poisson:        mu log
    "13" = 2L,   # hurdle_neg_binomial_2: mu log
    gdpar_abort(
      sprintf(
        "Internal error: no canonical inv_link id registered for stan_id %s.",
        as.character(stan_id)
      ),
      class = "gdpar_internal_error",
      data = list(stan_id = stan_id)
    )
  )
}

#' Internal: canonical inverse-link ID per built-in family for an
#' arbitrary slot index
#'
#' Returns the integer inverse-link ID of the slot with the supplied
#' index in the canonical param_specs of the family with the supplied
#' stan_id. Used by the homogeneous branch of the heterogeneous
#' resolver: when family_id_k[k] == family_id_k[1], the slot k uses
#' the canonical link of slot k of the location family (NOT the
#' canonical link of slot 1 via L1).
#'
#' @param stan_id Integer scalar.
#' @param slot_idx Integer scalar; 1-based slot index in the family's
#'   canonical param_specs.
#' @return Integer scalar in {0, 1, 2}.
#' @keywords internal
#' @noRd
.gdpar_canonical_inv_link_id_slot <- function(stan_id, slot_idx) {
  if (is.null(stan_id) || is.na(stan_id)) {
    gdpar_abort(
      "Internal error: canonical inv_link id requested for NA stan_id.",
      class = "gdpar_internal_error"
    )
  }
  # stan_id 7 (lognormal_loc_scale) is a custom_K registry pattern, not
  # a built-in gdpar_family. Its slot specs live in
  # .gdpar_K_custom_patterns(); look them up directly to avoid the
  # gdpar_family / gdpar_family_custom_K branch divergence.
  if (as.integer(stan_id) == 7L) {
    patterns <- .gdpar_K_custom_patterns()
    pattern <- patterns[["lognormal_loc_scale"]]
    if (slot_idx < 1L || slot_idx > length(pattern$slot_specs)) {
      gdpar_abort(
        sprintf(
          paste(
            "Internal error: slot_idx %d out of range for custom_K pattern",
            "'lognormal_loc_scale' (1:%d)."
          ),
          slot_idx, length(pattern$slot_specs)
        ),
        class = "gdpar_internal_error",
        data = list(stan_id = 7L, slot_idx = slot_idx,
                    n_slots = length(pattern$slot_specs))
      )
    }
    link <- pattern$slot_specs[[slot_idx]]$link
    return(switch(link,
      identity = 0L,
      logit    = 1L,
      log      = 2L,
      gdpar_abort(
        sprintf(
          "Internal error: link '%s' not supported by apply_inv_link_by_id.",
          link
        ),
        class = "gdpar_internal_error",
        data = list(link = link)
      )
    ))
  }
  family_name <- switch(as.character(as.integer(stan_id)),
    "1"  = "gaussian",
    "3"  = "neg_binomial_2",
    "5"  = "beta",
    "6"  = "gamma",
    "8"  = "student_t",
    "9"  = "tweedie",
    "10" = "zip",
    "11" = "zinb",
    "12" = "hurdle_poisson",
    "13" = "hurdle_neg_binomial_2",
    NULL
  )
  if (is.null(family_name)) {
    gdpar_abort(
      sprintf(
        "Internal error: no built-in family registered for stan_id %s.",
        as.character(stan_id)
      ),
      class = "gdpar_internal_error",
      data = list(stan_id = stan_id)
    )
  }
  specs <- .gdpar_family_param_specs_for(family_name,
                                          .gdpar_default_link_for(family_name))
  if (slot_idx < 1L || slot_idx > length(specs)) {
    gdpar_abort(
      sprintf(
        "Internal error: slot_idx %d out of range for family '%s' (1:%d).",
        slot_idx, family_name, length(specs)
      ),
      class = "gdpar_internal_error",
      data = list(stan_id = stan_id, slot_idx = slot_idx,
                  n_slots = length(specs))
    )
  }
  link <- specs[[slot_idx]]$link
  switch(link,
    identity = 0L,
    logit    = 1L,
    log      = 2L,
    gdpar_abort(
      sprintf(
        "Internal error: link '%s' not supported by apply_inv_link_by_id.",
        link
      ),
      class = "gdpar_internal_error",
      data = list(link = link)
    )
  )
}

#' Internal: default link for a built-in family by name
#'
#' Used by .gdpar_canonical_inv_link_id_slot to look up param_specs
#' without requiring the caller to know the link conventions.
#'
#' @keywords internal
#' @noRd
.gdpar_default_link_for <- function(name) {
  switch(name,
    gaussian              = "identity",
    poisson               = "log",
    neg_binomial_2        = "log",
    bernoulli             = "logit",
    beta                  = "logit",
    gamma                 = "log",
    student_t             = "identity",
    tweedie               = "log",
    zip                   = "log",
    zinb                  = "log",
    hurdle_poisson        = "log",
    hurdle_neg_binomial_2 = "log"
  )
}

#' Internal: canonical support of slot 1 (location) of a built-in family
#'
#' Returns the canonical support string emitted by the location slot
#' of the family with the supplied stan_id. Used by the heterogeneous
#' validator (D4) of 8.3.7: when family_id_k[k] is heterogeneous, the
#' emitted support is the support of slot 1 of that family per L1.
#'
#' @param stan_id Integer scalar.
#' @return Character scalar; one of "real_line", "positive_real",
#'   "unit_interval", "bounded_open".
#' @keywords internal
#' @noRd
.gdpar_canonical_support_slot1 <- function(stan_id) {
  switch(as.character(as.integer(stan_id)),
    "1"  = "real_line",
    "3"  = "positive_real",
    "5"  = "unit_interval",
    "6"  = "positive_real",
    "7"  = "real_line",
    "8"  = "real_line",
    "9"  = "positive_real",
    "10" = "positive_real",
    "11" = "positive_real",
    "12" = "positive_real",
    "13" = "positive_real",
    gdpar_abort(
      sprintf(
        "Internal error: no canonical support registered for stan_id %s.",
        as.character(stan_id)
      ),
      class = "gdpar_internal_error",
      data = list(stan_id = stan_id)
    )
  )
}

#' Internal: subset coherence predicate between supports
#'
#' Returns TRUE when the support emitted by a heterogeneous family
#' (slot 1) is a coherent subset of the support required by the role
#' of the slot in the location family. The relation is defined as:
#' \itemize{
#'   \item real_line is subset of \{real_line\}.
#'   \item positive_real is subset of \{positive_real, real_line\}.
#'   \item unit_interval is subset of \{unit_interval, positive_real,
#'     real_line\}.
#'   \item bounded_open is subset of \{bounded_open, positive_real,
#'     real_line\} (without checking bounds: 8.3.7 does not register
#'     bounded ranges for slot 1 of any built-in family).
#' }
#'
#' Used by the D4 validator of 8.3.7. The reflexive case
#' (emitted == required) always passes. Custom supports (\code{custom})
#' are rejected as the heterogeneous family cannot have a custom slot 1
#' built-in.
#'
#' @param emitted Character scalar with the support emitted by the
#'   heterogeneous family's slot 1.
#' @param required Character scalar with the support required by the
#'   role of the slot in the location family.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.gdpar_support_subset_coherent <- function(emitted, required) {
  if (identical(emitted, required)) {
    return(TRUE)
  }
  if (identical(emitted, "real_line")) {
    return(identical(required, "real_line"))
  }
  if (identical(emitted, "positive_real")) {
    return(required %in% c("positive_real", "real_line"))
  }
  if (identical(emitted, "unit_interval")) {
    return(required %in% c("unit_interval", "positive_real", "real_line"))
  }
  if (identical(emitted, "bounded_open")) {
    return(required %in% c("bounded_open", "positive_real", "real_line"))
  }
  FALSE
}

#' Internal: enumerate alternative built-in families with a compatible
#' slot 1 support
#'
#' Used by the D4 validator to build informative error messages:
#' when a user-supplied heterogeneous family fails the subset coherence
#' predicate, the message suggests built-in families whose slot 1
#' support would satisfy the requirement.
#'
#' Excludes the location family itself (homogeneous trivial case is
#' not a suggestion) and the family that already failed.
#'
#' @param required Character scalar; the support required by the role.
#' @param exclude_stan_ids Integer vector of stan_ids to exclude.
#' @return Character vector of family names (built-in only).
#' @keywords internal
#' @noRd
.gdpar_compatible_families_for_support <- function(required,
                                                    exclude_stan_ids = integer(0)) {
  candidates <- list(
    list(name = "gaussian",              stan_id = 1L),
    list(name = "neg_binomial_2",        stan_id = 3L),
    list(name = "beta",                  stan_id = 5L),
    list(name = "gamma",                 stan_id = 6L),
    list(name = "student_t",             stan_id = 8L),
    list(name = "tweedie",               stan_id = 9L)
  )
  out <- character(0)
  for (c in candidates) {
    if (c$stan_id %in% exclude_stan_ids) next
    emitted <- .gdpar_canonical_support_slot1(c$stan_id)
    if (.gdpar_support_subset_coherent(emitted, required)) {
      out <- c(out, c$name)
    }
  }
  out
}

#' Internal: validate a heterogeneous family list against the location
#' family's slot roles
#'
#' Implements D4 (validator strict + per param_spec$support + rich
#' errors + did_status from canonical primary) of the 8.3.7 scoping.
#' Iterates over the slots k = 2:K of the named-list family and
#' checks the subset coherence predicate between the emitted support
#' (slot 1 of family_het_list[[k]]) and the required support (slot k
#' of the location family's promoted param_specs).
#'
#' Slot 1 (location) is not validated here: it is the family that
#' governs the likelihood and determines the role assignment of the
#' remaining slots. The dispatcher in gdpar() must pass the location
#' family separately as the canonical reference.
#'
#' On the first violation, aborts with gdpar_input_error citing the
#' offending slot name, role, emitted support, required support, the
#' family name of the location family, and a suggestion list of
#' built-in families whose slot 1 support would satisfy the predicate.
#'
#' @param family_het_list Named list of gdpar_family objects, length K.
#'   Names must match slot_names. Entry [[1]] is the location family.
#' @param location_param_specs List of K param_specs from the location
#'   family after promotion to per_observation (the canonical role
#'   assignment of the K slots).
#' @param slot_names Character vector of length K with the canonical
#'   slot names.
#' @return Invisible NULL on success.
#' @keywords internal
#' @noRd
.gdpar_validate_heterogeneous_family_K <- function(family_het_list,
                                                    location_param_specs,
                                                    slot_names) {
  K <- length(family_het_list)
  if (K != length(location_param_specs) || K != length(slot_names)) {
    gdpar_abort(
      sprintf(
        paste(
          "Internal error: heterogeneous validator received inconsistent K:",
          "family_het_list = %d, location_param_specs = %d, slot_names = %d."
        ),
        K, length(location_param_specs), length(slot_names)
      ),
      class = "gdpar_internal_error"
    )
  }
  location_family <- family_het_list[[1L]]
  location_stan_id <- as.integer(location_family$stan_id)
  for (k in seq_len(K)[-1L]) {
    f_k <- family_het_list[[k]]
    if (!inherits(f_k, "gdpar_family")) {
      gdpar_abort(
        sprintf(
          paste(
            "Heterogeneous family entry for slot '%s' must be a gdpar_family",
            "object; received an object of class %s."
          ),
          slot_names[k],
          paste(sQuote(class(f_k)), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_names[k],
                    received_class = class(f_k))
      )
    }
    f_k_stan_id <- f_k$stan_id
    if (is.null(f_k_stan_id) || is.na(f_k_stan_id)) {
      gdpar_abort(
        sprintf(
          paste(
            "Heterogeneous family entry for slot '%s' has stan_id NA (custom",
            "free-form family). Sub-phase 8.3.7 admits only built-in families",
            "and descriptor-based custom_K patterns as heterogeneous slot",
            "entries."
          ),
          slot_names[k]
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_names[k], family = f_k$name)
      )
    }
    f_k_stan_id <- as.integer(f_k_stan_id)
    if (f_k_stan_id == location_stan_id) {
      next
    }
    emitted <- .gdpar_canonical_support_slot1(f_k_stan_id)
    required <- location_param_specs[[k]]$support
    if (.gdpar_support_subset_coherent(emitted, required)) {
      next
    }
    role <- location_param_specs[[k]]$family_role
    suggestions <- .gdpar_compatible_families_for_support(
      required = required,
      exclude_stan_ids = c(location_stan_id, f_k_stan_id)
    )
    sugg_msg <- if (length(suggestions) > 0L) {
      paste(sQuote(suggestions), collapse = ", ")
    } else {
      "(no built-in family has a slot 1 support compatible with this role)"
    }
    gdpar_abort(
      sprintf(
        paste(
          "Heterogeneous family for slot '%s' (role '%s' in location",
          "family '%s') emits canonical support '%s', which is not a",
          "subset of the required support '%s'. Suggested alternatives:",
          "%s."
        ),
        slot_names[k], role, location_family$name,
        emitted, required, sugg_msg
      ),
      class = "gdpar_input_error",
      data = list(slot          = slot_names[k],
                  role          = role,
                  emitted       = emitted,
                  required      = required,
                  location_family = location_family$name,
                  het_family    = f_k$name,
                  suggestions   = suggestions)
    )
  }
  invisible(NULL)
}

#' Internal: compute the inv_link_id_per_slot vector for a fit
#'
#' Returns an integer vector of length K with the inverse-link ID per
#' slot under the L1 refined rule (D3.5 of 8.3.7):
#' \itemize{
#'   \item If family_id_k[k] == family_id_k[1] (homogeneous slot k),
#'     return the canonical inv_link ID of slot k of the location
#'     family.
#'   \item Otherwise (heterogeneous slot k), return the canonical
#'     inv_link ID of slot 1 of the family of slot k (L1).
#' }
#'
#' The vector is passed to the Stan template via the data field
#' inv_link_id_per_slot. The helper apply_inv_link_by_id() in the
#' template dispatches on this ID at the model and generated quantities
#' blocks.
#'
#' Backward-compat note (D7 of 8.3.7): in the strict homogeneous
#' regime (all family_id_k entries equal), the vector reproduces the
#' canonical slot-by-slot links of the location family. The Stan
#' refactor of the likelihood branches uses apply_inv_link_by_id() in
#' place of the previously hardcoded inv_link calls, which preserves
#' mathematical equivalence but is NOT guaranteed bit-for-bit
#' identical with pre-8.3.7 outputs. Goldens K=2 Beta and Gamma must
#' be re-bootstrapped on close.
#'
#' @param family_id_k_vector Integer vector of length K with the
#'   per-slot stan_id.
#' @param location_family The location family (slot 1), used to look
#'   up the canonical slot-k links in the homogeneous case.
#' @return Integer vector of length K with values in {0, 1, 2}.
#' @keywords internal
#' @noRd
.gdpar_compute_inv_link_id_per_slot <- function(family_id_k_vector,
                                                 location_family) {
  K <- length(family_id_k_vector)
  out <- integer(K)
  location_stan_id <- as.integer(location_family$stan_id)
  for (k in seq_len(K)) {
    f_k_id <- as.integer(family_id_k_vector[k])
    if (f_k_id == location_stan_id) {
      out[k] <- .gdpar_canonical_inv_link_id_slot(location_stan_id, k)
    } else {
      out[k] <- .gdpar_canonical_inv_link_id_slot1(f_k_id)
    }
  }
  out
}

#' Internal: resolve a heterogeneous named-list family argument
#'
#' Implements the named-list public API (D5 of 8.3.7) and the
#' materialization rule (D3.5 / L1 refined). Returns a list with
#' \code{location_family} (the gdpar_family that governs the
#' likelihood, with per-slot did_status overrides propagated from the
#' heterogeneous slot families per L1), \code{family_id_k_vector}
#' (length-K integer vector with the stan_id of each slot), and
#' \code{is_heterogeneous} (TRUE if at least one slot has a different
#' stan_id from slot 1).
#'
#' The location_family returned has its param_specs already updated:
#' for each heterogeneous slot k (k >= 2 with family_id_k[k] !=
#' family_id_k[1]), the slot's did_status, did_condition and
#' did_reference fields are replaced with those of slot 1 of the
#' heterogeneous family, per L1. The link, support and
#' prior_canonical_kind of the slot are NOT overwritten in
#' param_specs (those fields stay as the original canonical role
#' assignment of the location family) because the Stan side computes
#' the effective inv_link via inv_link_id_per_slot and the priors are
#' applied via the canonical THETA_REF_PRIOR_BLOCK on the
#' linear-predictor space, which is invariant under the link choice
#' of the slot. (Re-evaluation of canonical prior overrides per slot
#' is deferred to Session 8.4 per gotcha 55.)
#'
#' @param family_input Either a gdpar_family object (homogeneous,
#'   passed-through) or a named list of gdpar_family objects keyed by
#'   slot name (heterogeneous).
#' @param slot_names Character vector of length K with the canonical
#'   slot names of the gdpar_formula_set.
#' @return Named list with components \code{location_family},
#'   \code{family_id_k_vector}, \code{is_heterogeneous}.
#' @keywords internal
#' @noRd
.gdpar_resolve_heterogeneous_family_K <- function(family_input, slot_names) {
  K <- length(slot_names)
  if (inherits(family_input, "gdpar_family")) {
    family_id_k_vector <- rep(as.integer(family_input$stan_id), K)
    return(list(
      location_family    = family_input,
      family_id_k_vector = family_id_k_vector,
      is_heterogeneous   = FALSE
    ))
  }
  if (!is.list(family_input)) {
    gdpar_abort(
      paste(
        "Argument 'family' must be a gdpar_family object or a named list of",
        "gdpar_family objects (one per slot) for the heterogeneous family",
        "path of sub-phase 8.3.7."
      ),
      class = "gdpar_input_error",
      data = list(received_class = class(family_input))
    )
  }
  if (is.null(names(family_input)) ||
      any(!nzchar(names(family_input))) ||
      anyDuplicated(names(family_input)) > 0L) {
    gdpar_abort(
      paste(
        "When 'family' is a named list (heterogeneous path), every entry",
        "must carry a unique non-empty slot name matching the slot names",
        "of the gdpar_formula_set."
      ),
      class = "gdpar_input_error"
    )
  }
  if (length(family_input) != K) {
    gdpar_abort(
      sprintf(
        paste(
          "Named-list family argument has %d entries but the formula set",
          "declares K = %d slots (%s)."
        ),
        length(family_input), K,
        paste(sQuote(slot_names), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(K = K, received_K = length(family_input),
                  slot_names = slot_names)
    )
  }
  missing_slots <- setdiff(slot_names, names(family_input))
  extra_slots   <- setdiff(names(family_input), slot_names)
  if (length(missing_slots) > 0L || length(extra_slots) > 0L) {
    gdpar_abort(
      sprintf(
        paste(
          "Named-list family slot names do not match the formula set.",
          "Missing slot(s): %s. Unexpected slot name(s): %s. Required",
          "slot order: %s."
        ),
        if (length(missing_slots) > 0L)
          paste(sQuote(missing_slots), collapse = ", ") else "(none)",
        if (length(extra_slots) > 0L)
          paste(sQuote(extra_slots), collapse = ", ") else "(none)",
        paste(sQuote(slot_names), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(missing = missing_slots, extra = extra_slots,
                  required = slot_names)
    )
  }
  family_het_list <- family_input[slot_names]
  for (k in seq_len(K)) {
    if (!inherits(family_het_list[[k]], "gdpar_family")) {
      gdpar_abort(
        sprintf(
          paste(
            "Entry '%s' of the named-list family must be a gdpar_family",
            "object; received an object of class %s."
          ),
          slot_names[k],
          paste(sQuote(class(family_het_list[[k]])), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_names[k],
                    received_class = class(family_het_list[[k]]))
      )
    }
  }
  location_family <- family_het_list[[1L]]
  if (is.null(location_family$stan_id) || is.na(location_family$stan_id)) {
    gdpar_abort(
      sprintf(
        paste(
          "The location-slot family '%s' has stan_id NA (custom free-form",
          "family). Sub-phase 8.3.7 requires the location slot family to",
          "be a built-in or descriptor-based custom_K family."
        ),
        location_family$name
      ),
      class = "gdpar_input_error",
      data = list(family = location_family$name)
    )
  }
  location_promoted <- .gdpar_promote_scope_per_observation(
    location_family, slot_names
  )
  location_stan_id <- as.integer(location_promoted$stan_id)
  family_id_k_vector <- unname(vapply(family_het_list,
                                function(f) as.integer(f$stan_id),
                                integer(1L)))
  is_heterogeneous <- any(family_id_k_vector != location_stan_id)
  heterogeneity_supported_locations <- c(1L, 3L, 5L, 6L, 7L)
  if (is_heterogeneous &&
      !(location_stan_id %in% heterogeneity_supported_locations)) {
    location_name_map <- c(`1` = "gaussian", `3` = "neg_binomial_2",
                            `5` = "beta", `6` = "gamma",
                            `7` = "lognormal_loc_scale (custom_K)",
                            `8` = "student_t", `9` = "tweedie",
                            `10` = "zip", `11` = "zinb",
                            `12` = "hurdle_poisson",
                            `13` = "hurdle_neg_binomial_2")
    supported_names <- location_name_map[
      as.character(heterogeneity_supported_locations)
    ]
    gdpar_abort(
      sprintf(
        paste(
          "Heterogeneous families per slot are wired in sub-phase 8.3.7 only",
          "for the K = 2 location families %s. The supplied location family",
          "'%s' (stan_id = %d) is a K = 3+ family whose Stan-side dispatch",
          "retains specialised idioms (*_log_lpmf, bernoulli_logit_lpmf,",
          "custom tweedie_lpdf); heterogeneous slot k >= 2 for this family",
          "is queued for a future sub-phase. Pass a single gdpar_family",
          "object instead of a named list to keep the homogeneous regime."
        ),
        paste(sQuote(supported_names), collapse = ", "),
        location_family$name, location_stan_id
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(location_family = location_family$name,
                  location_stan_id = location_stan_id,
                  supported_locations = heterogeneity_supported_locations)
    )
  }
  .gdpar_validate_heterogeneous_family_K(
    family_het_list      = family_het_list,
    location_param_specs = location_promoted$param_specs,
    slot_names           = slot_names
  )
  for (k in seq_len(K)[-1L]) {
    f_k_id <- family_id_k_vector[k]
    if (f_k_id == location_stan_id) next
    f_k <- family_het_list[[k]]
    location_promoted$param_specs[[k]]$did_status    <- f_k$did_status
    location_promoted$param_specs[[k]]$did_condition <- f_k$did_condition
    location_promoted$param_specs[[k]]$did_reference <- f_k$did_reference
  }
  list(
    location_family    = location_promoted,
    family_id_k_vector = family_id_k_vector,
    is_heterogeneous   = is_heterogeneous
  )
}

#' Internal: detect intercept suppression in a formula
#'
#' Detects the patterns `- 1` and `+ 0` in either side of a formula by
#' inspecting the intercept attribute of \code{\link[stats]{terms}}. In
#' the AMM canonical form the population anchor \eqn{\theta_{ref}} is a
#' structural term and cannot be suppressed: a user request to remove it
#' via R's intercept-suppression syntax is interpreted as an abdication
#' of the canonical contract and aborts with \code{gdpar_input_error}.
#' Decision recorded in
#' \code{project_gdpar_block_8_3_extended_plan} and
#' \code{project_gdpar_block_8_session_1_decisions} (memory entries,
#' 2026-05-20).
#'
#' @param f A formula object (one-sided or two-sided).
#' @param slot_name Character scalar with the slot name from the
#'   enclosing \code{gdpar_formula_set} used in the error message.
#'
#' @return Invisibly TRUE when the formula keeps its intercept;
#'   otherwise aborts.
#' @keywords internal
#' @noRd
.gdpar_check_no_intercept_suppression <- function(f, slot_name) {
  trm <- tryCatch(stats::terms(f),
                  error = function(e) NULL)
  if (is.null(trm)) {
    invisible(TRUE)
    return(invisible(TRUE))
  }
  if (isTRUE(attr(trm, "intercept") == 0L)) {
    gdpar_abort(
      sprintf(
        paste(
          "Slot '%s' of the formula set contains intercept suppression",
          "(`- 1` or `+ 0`): %s. AMM canoniza el anchor poblacional",
          "theta_ref como termino estructural; suprimirlo via `- 1` o",
          "`+ 0` esta fuera del modelo. Para modelos sin anchor",
          "poblacional, usar glm()/brms/stan_glm() directamente."
        ),
        slot_name, deparse(f, width.cutoff = 500L)
      ),
      class = "gdpar_input_error",
      data = list(slot = slot_name,
                  formula = deparse(f, width.cutoff = 500L))
    )
  }
  invisible(TRUE)
}

#' Internal: assert that all arguments to a constructor are named formulas
#'
#' @param args A list captured from \code{...}.
#' @param fn_name Character scalar with the calling function's name for
#'   inclusion in the error message.
#' @return Invisibly TRUE on success.
#' @keywords internal
#' @noRd
.gdpar_assert_named_formula_list <- function(args, fn_name) {
  if (length(args) == 0L) {
    gdpar_abort(
      sprintf("%s() requires at least one named formula.", fn_name),
      class = "gdpar_input_error"
    )
  }
  arg_names <- names(args)
  if (is.null(arg_names) || any(!nzchar(arg_names))) {
    gdpar_abort(
      sprintf(
        paste(
          "All formulas passed to %s() must be named with the canonical",
          "parameter name (e.g., mu = y ~ a(x), sigma = ~ a(z))."
        ),
        fn_name
      ),
      class = "gdpar_input_error"
    )
  }
  if (anyDuplicated(arg_names) > 0L) {
    dups <- unique(arg_names[duplicated(arg_names)])
    gdpar_abort(
      sprintf(
        "Duplicate slot names in %s(): %s.",
        fn_name, paste(sQuote(dups), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(duplicates = dups)
    )
  }
  for (i in seq_along(args)) {
    if (!inherits(args[[i]], "formula")) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' of %s() must be a formula; received an object",
            "of class %s."
          ),
          arg_names[[i]], fn_name,
          paste(sQuote(class(args[[i]])), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(slot = arg_names[[i]],
                    received = class(args[[i]]))
      )
    }
  }
  invisible(TRUE)
}

#' Construct a canonical formula set for multi-parameter AMM modeling
#'
#' Build the canonical internal representation that \code{\link{gdpar}}
#' consumes when more than one structural parameter of the family is
#' modeled with an AMM design (multi-parametric distributional regression
#' in the sense of decision D-F1 of Block 8 Session 1, materialized in
#' sub-phase 8.3.3 of the package). One slot per individual parameter,
#' named with the canonical parameter name as declared by the family's
#' \code{param_specs}; the first slot carries the outcome variable on
#' its left-hand side, while subsequent slots are one-sided formulas
#' that describe the AMM design of an auxiliary parameter (e.g., the
#' Gaussian \code{sigma} or the negative-binomial \code{phi}).
#'
#' This constructor is the canonical low-cost entry point; the
#' brms-style sugar \code{\link{gdpar_bf}} produces an equivalent object
#' from a sequence of two-sided formulas.
#'
#' @section Validation contract:
#' \itemize{
#'   \item All formulas must be named with the canonical parameter
#'     name (no positional formulas).
#'   \item The first slot must be a two-sided formula whose left-hand
#'     side is a single symbol naming the outcome variable
#'     (e.g., \code{mu = y ~ a(x)}). Function calls on the LHS such
#'     as \code{log(y)} are rejected to keep the contract explicit;
#'     pre-transform the outcome in the data frame instead.
#'   \item Subsequent slots must be one-sided formulas
#'     (e.g., \code{sigma = ~ a(x)}).
#'   \item Slot names must be unique.
#'   \item No formula may suppress the intercept via \code{-1} or
#'     \code{+0}: the AMM population anchor \eqn{\theta_{ref}} is a
#'     structural term that cannot be removed. The constructor aborts
#'     with an informative \code{gdpar_input_error} when this pattern
#'     is detected.
#' }
#'
#' Set membership against the family's eligible \code{param_specs}
#' (slot names must be a subset of the eligible parameters of the
#' family passed to \code{\link{gdpar}}) is enforced downstream by
#' \code{\link{gdpar}} once both the formula set and the family are
#' known.
#'
#' @section Methodological notes:
#' Block 8 Session 1 decision 1C established that a family is a list
#' of parameter specifications. Block 8.3.3 closes the high-level API:
#' the formula set is the single source of truth for the K-individual
#' parameter set; any parameter named in this object is promoted to
#' \code{scope = "per_observation"} in the family copy that
#' \code{\link{gdpar}} uses internally. Parameters of the family that
#' are not named in the formula set remain with their canonical
#' \code{scope} (typically \code{population}). See memory entries
#' \code{project_gdpar_block_8_3_extended_plan} and
#' \code{project_gdpar_block_8_session_1_decisions} (2026-05-20).
#'
#' @param ... Named formulas. The first must be two-sided
#'   (\code{mu = y ~ a(x)}). Subsequent ones must be one-sided
#'   (\code{sigma = ~ a(x)}).
#'
#' @return An object of class \code{gdpar_formula_set} with components
#'   \code{outcome} (character scalar with the outcome variable name),
#'   \code{formulas} (named list of formula objects; the first is
#'   two-sided, subsequent ones are one-sided), \code{param_names}
#'   (character vector of slot names, identical to
#'   \code{names(formulas)}) and \code{env} (the environment of the
#'   first formula, used downstream for evaluation).
#'
#' @section Dependencies:
#' Uses \code{\link[stats]{terms}} to detect intercept suppression.
#'
#' @examples
#' fs <- gdpar_formula_set(
#'   mu    = y ~ a(x1) + b(z1),
#'   sigma = ~ a(x2)
#' )
#' print(fs)
#' fs[["mu"]]
#' names(fs)
#'
#' @seealso \code{\link{gdpar_bf}}, \code{\link{gdpar}},
#'   \code{\link{gdpar_family}}
#' @export
gdpar_formula_set <- function(...) {
  args <- list(...)
  .gdpar_assert_named_formula_list(args, "gdpar_formula_set")
  arg_names <- names(args)
  first_f <- args[[1L]]
  if (length(first_f) != 3L) {
    gdpar_abort(
      sprintf(
        paste(
          "The first slot ('%s') of gdpar_formula_set() must be a",
          "two-sided formula whose left-hand side is the outcome",
          "variable (e.g., %s = y ~ a(x)). Received: %s."
        ),
        arg_names[[1L]], arg_names[[1L]],
        deparse(first_f, width.cutoff = 500L)
      ),
      class = "gdpar_input_error",
      data = list(slot = arg_names[[1L]],
                  received = deparse(first_f, width.cutoff = 500L))
    )
  }
  outcome_expr <- first_f[[2L]]
  if (!is.symbol(outcome_expr)) {
    gdpar_abort(
      sprintf(
        paste(
          "The left-hand side of the first slot ('%s') of",
          "gdpar_formula_set() must be a single variable name; received: %s.",
          "Pre-transform the outcome in the data frame if a function of y",
          "is needed."
        ),
        arg_names[[1L]], deparse(outcome_expr, width.cutoff = 500L)
      ),
      class = "gdpar_input_error",
      data = list(slot = arg_names[[1L]],
                  received = deparse(outcome_expr, width.cutoff = 500L))
    )
  }
  outcome <- as.character(outcome_expr)
  if (length(args) > 1L) {
    for (i in 2:length(args)) {
      f <- args[[i]]
      if (length(f) != 2L) {
        gdpar_abort(
          sprintf(
            paste(
              "Slot '%s' of gdpar_formula_set() must be a one-sided formula",
              "(e.g., %s = ~ a(x)). Two-sided formulas are accepted only in",
              "the first slot, which carries the outcome on its LHS.",
              "Received: %s."
            ),
            arg_names[[i]], arg_names[[i]],
            deparse(f, width.cutoff = 500L)
          ),
          class = "gdpar_input_error",
          data = list(slot = arg_names[[i]],
                      received = deparse(f, width.cutoff = 500L))
        )
      }
    }
  }
  for (i in seq_along(args)) {
    .gdpar_check_no_intercept_suppression(args[[i]], arg_names[[i]])
  }
  obj <- list(
    outcome     = outcome,
    formulas    = args,
    param_names = arg_names,
    env         = environment(first_f)
  )
  class(obj) <- c("gdpar_formula_set", "list")
  obj
}

#' Construct a formula set with brms-style sugar
#'
#' Build a \code{\link{gdpar_formula_set}} from a sequence of two-sided
#' formulas in the style of \code{brms::bf}. The first formula carries
#' the outcome on its LHS and defaults the first slot name to
#' \code{"mu"} (the canonical name of the location parameter across
#' built-in families); each subsequent formula must carry the canonical
#' parameter name as its LHS (e.g., \code{sigma ~ a(x)}). The result is
#' an object of class \code{gdpar_formula_set} identical to what the
#' explicit constructor produces.
#'
#' This is sugar: \code{gdpar_bf(y ~ a(x1), sigma ~ a(x2))} is
#' equivalent to
#' \code{gdpar_formula_set(mu = y ~ a(x1), sigma = ~ a(x2))}.
#' Naming the first argument overrides the default \code{"mu"}:
#' \code{gdpar_bf(theta = y ~ a(x))} produces a single-slot set with
#' name \code{"theta"}, intended for custom families whose location
#' parameter is not canonically called \code{mu}.
#'
#' The brms package is not a runtime dependency of \code{gdpar}; the
#' constructor returns a native \code{gdpar_formula_set} regardless of
#' whether \pkg{brms} is installed.
#'
#' @param ... Two-sided formulas. The first carries the outcome on its
#'   LHS; subsequent ones carry the canonical parameter name on their
#'   LHS.
#'
#' @return An object of class \code{gdpar_formula_set} (see
#'   \code{\link{gdpar_formula_set}} for components).
#'
#' @section Dependencies:
#' Internally calls \code{\link{gdpar_formula_set}}.
#'
#' @examples
#' fs <- gdpar_bf(y ~ a(x1) + b(z1), sigma ~ a(x2))
#' print(fs)
#'
#' @seealso \code{\link{gdpar_formula_set}}, \code{\link{gdpar}}
#' @export
gdpar_bf <- function(...) {
  args <- list(...)
  if (length(args) == 0L) {
    gdpar_abort(
      "gdpar_bf() requires at least one formula.",
      class = "gdpar_input_error"
    )
  }
  for (i in seq_along(args)) {
    if (!inherits(args[[i]], "formula")) {
      gdpar_abort(
        sprintf(
          paste(
            "Argument %d of gdpar_bf() must be a formula; received an",
            "object of class %s."
          ),
          i, paste(sQuote(class(args[[i]])), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(position = i, received = class(args[[i]]))
      )
    }
  }
  arg_names <- names(args)
  if (is.null(arg_names)) arg_names <- rep("", length(args))
  first_f <- args[[1L]]
  if (length(first_f) != 3L) {
    gdpar_abort(
      sprintf(
        paste(
          "The first formula of gdpar_bf() must be two-sided with the",
          "outcome variable on its left-hand side (e.g., y ~ a(x)).",
          "Received: %s."
        ),
        deparse(first_f, width.cutoff = 500L)
      ),
      class = "gdpar_input_error",
      data = list(received = deparse(first_f, width.cutoff = 500L))
    )
  }
  slot_names <- character(length(args))
  slot_names[[1L]] <- if (nzchar(arg_names[[1L]])) arg_names[[1L]] else "mu"
  if (length(args) > 1L) {
    for (i in 2:length(args)) {
      f <- args[[i]]
      if (length(f) != 3L) {
        gdpar_abort(
          sprintf(
            paste(
              "Formula %d of gdpar_bf() must be two-sided with the parameter",
              "name on its left-hand side (e.g., sigma ~ a(x)). Received: %s."
            ),
            i, deparse(f, width.cutoff = 500L)
          ),
          class = "gdpar_input_error",
          data = list(position = i,
                      received = deparse(f, width.cutoff = 500L))
        )
      }
      lhs <- f[[2L]]
      if (!is.symbol(lhs)) {
        gdpar_abort(
          sprintf(
            paste(
              "The left-hand side of formula %d of gdpar_bf() must be a",
              "single parameter name; received: %s."
            ),
            i, deparse(lhs, width.cutoff = 500L)
          ),
          class = "gdpar_input_error",
          data = list(position = i,
                      received = deparse(lhs, width.cutoff = 500L))
        )
      }
      explicit_name <- if (nzchar(arg_names[[i]])) arg_names[[i]] else NA_character_
      lhs_name <- as.character(lhs)
      if (!is.na(explicit_name) && explicit_name != lhs_name) {
        gdpar_abort(
          sprintf(
            paste(
              "Formula %d of gdpar_bf() has an explicit argument name '%s'",
              "that does not match its LHS '%s'. Pass either the named form",
              "with no LHS (e.g., %s = ~ a(x)) or the positional form with",
              "the parameter name on the LHS (e.g., %s ~ a(x))."
            ),
            i, explicit_name, lhs_name, explicit_name, lhs_name
          ),
          class = "gdpar_input_error",
          data = list(position = i,
                      arg_name = explicit_name, lhs = lhs_name)
        )
      }
      slot_names[[i]] <- lhs_name
      one_sided <- stats::as.formula(
        paste("~", deparse(f[[3L]], width.cutoff = 500L)),
        env = environment(f)
      )
      args[[i]] <- one_sided
    }
  }
  if (anyDuplicated(slot_names) > 0L) {
    dups <- unique(slot_names[duplicated(slot_names)])
    gdpar_abort(
      sprintf(
        "Duplicate parameter names in gdpar_bf(): %s.",
        paste(sQuote(dups), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(duplicates = dups)
    )
  }
  names(args) <- slot_names
  do.call(gdpar_formula_set, args)
}

#' Print method for gdpar_formula_set objects
#'
#' @param x An object of class \code{gdpar_formula_set}.
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_formula_set <- function(x, ...) {
  K <- length(x$formulas)
  cat("<gdpar_formula_set> K = ", K, "\n", sep = "")
  cat("  outcome     : ", x$outcome, "\n", sep = "")
  cat("  param_names : ", paste(x$param_names, collapse = ", "),
      "\n", sep = "")
  cat("  formulas    :\n")
  for (k in x$param_names) {
    f <- x$formulas[[k]]
    cat("    ", k, " : ",
        deparse(f, width.cutoff = 500L), "\n", sep = "")
  }
  invisible(x)
}

#' Extract a single formula from a gdpar_formula_set
#'
#' @param x An object of class \code{gdpar_formula_set}.
#' @param i A character scalar with the slot name or an integer index.
#' @return The formula stored at slot \code{i}.
#' @export
`[[.gdpar_formula_set` <- function(x, i) {
  x$formulas[[i]]
}

#' Subset a gdpar_formula_set by slot name or position
#'
#' @param x An object of class \code{gdpar_formula_set}.
#' @param i A character vector of slot names or an integer vector of
#'   positions.
#' @return A named list of formulas (not a \code{gdpar_formula_set};
#'   downstream re-validation requires reconstruction via
#'   \code{\link{gdpar_formula_set}}).
#' @export
`[.gdpar_formula_set` <- function(x, i) {
  x$formulas[i]
}

#' Slot names of a gdpar_formula_set
#'
#' @param x An object of class \code{gdpar_formula_set}.
#' @return Character vector with the slot names in declaration order.
#' @export
names.gdpar_formula_set <- function(x) {
  x$param_names
}

#' Number of K-individual parameters in a gdpar_formula_set
#'
#' @param x An object of class \code{gdpar_formula_set}.
#' @return Integer scalar equal to the number of slots.
#' @export
length.gdpar_formula_set <- function(x) {
  length(x$formulas)
}

#' Internal: detect whether the RHS of a classic formula contains AMM wrappers
#'
#' Walks the right-hand side of a classic two-sided formula and returns
#' \code{TRUE} when at least one top-level summand is a call to one of
#' the AMM wrappers (\code{a}, \code{b}, \code{W}). The function is
#' non-throwing: an invalid RHS that the parser would later reject
#' still returns \code{TRUE} here when an AMM wrapper is detected, so
#' that the dispatch in \code{\link{gdpar}} routes the call to the
#' parser (which produces the canonical error message).
#'
#' Used by the dispatch decision P-dispatch in sub-phase 8.3.3:
#' classic formula + AMM wrappers in the RHS routes through the parser
#' (new K-individual path); classic formula without AMM wrappers routes
#' through the legacy single-\code{amm_spec} path bit-exactly.
#'
#' @param formula A classic two-sided formula.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.gdpar_rhs_has_amm_calls <- function(formula) {
  if (!inherits(formula, "formula")) return(FALSE)
  rhs <- if (length(formula) == 3L) formula[[3L]] else formula[[2L]]
  has_amm <- function(expr) {
    if (is.call(expr) && length(expr) == 3L &&
        identical(expr[[1L]], as.symbol("+"))) {
      return(has_amm(expr[[2L]]) || has_amm(expr[[3L]]))
    }
    if (is.call(expr) && is.symbol(expr[[1L]])) {
      nm <- as.character(expr[[1L]])
      if (nm %in% c("a", "b", "W")) return(TRUE)
    }
    FALSE
  }
  has_amm(rhs)
}

#' Internal: test whether an amm_spec object is the default Level 0 spec
#'
#' Returns \code{TRUE} when the object is an \code{amm_spec} with no
#' active components (\code{a}, \code{b}, \code{W} all \code{NULL},
#' no \code{x_vars}, \code{p == 1L} and no \code{dims}). Used by the
#' dispatch in \code{\link{gdpar}} to detect the canonical
#' \dQuote{user did not customize amm} case and treat it as an absent
#' argument when the new K-individual paths are taken.
#'
#' @param amm Any R object.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.gdpar_is_default_amm_spec <- function(amm) {
  inherits(amm, "amm_spec") &&
    is.null(amm$a) && is.null(amm$b) && is.null(amm$W) &&
    is.null(amm$x_vars) &&
    (is.null(amm$p) || isTRUE(amm$p == 1L)) &&
    is.null(amm$dims)
}

#' Internal: parse an AMM formula RHS into a sub-amm_spec descriptor
#'
#' Walks the right-hand side expression of a one-sided formula and
#' extracts the AMM components declared via the special function calls
#' \code{a(...)}, \code{b(...)} and \code{W()}. The output is a plain
#' list with components \code{a_formula} (one-sided formula or NULL),
#' \code{b_formula} (one-sided formula or NULL) and \code{W_present}
#' (logical scalar). Downstream the dispatch in \code{\link{gdpar}}
#' combines this list with the external \code{W_basis} argument to
#' build an \code{\link{amm_spec}} per slot.
#'
#' Scoping decisions canonized in memory entry
#' \code{project_gdpar_block_8_3_extended_plan} (2026-05-20):
#' \itemize{
#'   \item P1: unknown functions in the RHS abort with
#'     \code{gdpar_input_error} listing the permitted ones
#'     (\code{a}, \code{b}, \code{W}).
#'   \item P2: each of \code{a}, \code{b}, \code{W} may appear at most
#'     once per slot; multiple occurrences abort.
#'   \item P3: \code{W()} must be called with zero arguments; the
#'     modulating basis travels via the external \code{W} argument of
#'     \code{\link{gdpar}}.
#'   \item P4: the interior of \code{a(...)} and \code{b(...)} is an
#'     RHS-style expression that is wrapped in \code{~ <expr>} and
#'     passed downstream to \code{\link[stats]{model.matrix}}.
#'   \item P5: a bare RHS of \code{1} (intercept only) is accepted as
#'     a Level 0 AMM (anchor only, no additive/multiplicative/modulating
#'     component); any other bare RHS without AMM function calls aborts
#'     with an informative error.
#' }
#'
#' @param rhs_formula A one-sided formula whose RHS is to be parsed.
#' @param slot_name Character scalar with the slot name from the
#'   enclosing \code{gdpar_formula_set}, used for error messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{a_formula}}{A one-sided formula \code{~ <expr>} when
#'     \code{a(<expr>)} appears in the RHS; \code{NULL} otherwise.}
#'   \item{\code{b_formula}}{A one-sided formula \code{~ <expr>} when
#'     \code{b(<expr>)} appears in the RHS; \code{NULL} otherwise.}
#'   \item{\code{W_present}}{Logical scalar; \code{TRUE} when \code{W()}
#'     appears in the RHS.}
#' }
#'
#' @keywords internal
#' @noRd
.gdpar_parse_amm_formula <- function(rhs_formula, slot_name) {
  if (!inherits(rhs_formula, "formula")) {
    gdpar_abort(
      sprintf(
        "Internal error: .gdpar_parse_amm_formula expected a formula for slot '%s'.",
        slot_name
      ),
      class = "gdpar_internal_error"
    )
  }
  rhs_expr <- if (length(rhs_formula) == 3L) {
    rhs_formula[[3L]]
  } else {
    rhs_formula[[2L]]
  }
  env <- environment(rhs_formula)
  summands <- .gdpar_split_amm_summands(rhs_expr, slot_name)
  out <- list(a_formula = NULL, b_formula = NULL, W_present = FALSE)
  for (s in summands) {
    if (is.symbol(s) && identical(s, as.symbol("1"))) {
      next
    }
    if (is.numeric(s) && length(s) == 1L && s == 1) {
      next
    }
    if (!is.call(s)) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' of the formula set contains the bare term '%s' on its",
            "right-hand side. AMM canoniza el contrato: cada termino debe",
            "envolverse en uno de los wrappers `a()`, `b()` o `W()`; para una",
            "AMM degenerada (Level 0, anchor-only) usar `~ 1`."
          ),
          slot_name, deparse(s, width.cutoff = 500L)
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_name,
                    term = deparse(s, width.cutoff = 500L))
      )
    }
    head <- s[[1L]]
    if (!is.symbol(head)) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' contains a call whose head is not a simple function name: %s.",
            "Only the wrappers `a()`, `b()` and `W()` are admitted."
          ),
          slot_name, deparse(s, width.cutoff = 500L)
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_name,
                    term = deparse(s, width.cutoff = 500L))
      )
    }
    fname <- as.character(head)
    if (!fname %in% c("a", "b", "W")) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' contains the unknown function call `%s()` on its",
            "right-hand side. AMM canoniza los wrappers permitidos: `a()`",
            "(componente aditivo), `b()` (componente multiplicativo) y",
            "`W()` (componente modulante)."
          ),
          slot_name, fname
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_name, fname = fname,
                    allowed = c("a", "b", "W"))
      )
    }
    if (fname == "W") {
      if (length(s) > 1L) {
        gdpar_abort(
          sprintf(
            paste(
              "Slot '%s' calls `W(...)` with arguments. Per the canonization",
              "of sub-phase 8.3.3, `W()` must be invoked with zero",
              "arguments; the modulating basis travels via the external",
              "`W` argument of gdpar() (e.g., `W = W_basis(...)`)."
            ),
            slot_name
          ),
          class = "gdpar_input_error",
          data = list(slot = slot_name,
                      term = deparse(s, width.cutoff = 500L))
        )
      }
      if (isTRUE(out$W_present)) {
        gdpar_abort(
          sprintf(
            paste(
              "Slot '%s' calls `W()` more than once. Only one occurrence of",
              "each AMM wrapper is admitted per slot."
            ),
            slot_name
          ),
          class = "gdpar_input_error",
          data = list(slot = slot_name, fname = "W")
        )
      }
      out$W_present <- TRUE
      next
    }
    if (length(s) != 2L) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' calls `%s(...)` with %d arguments; it must be called",
            "with exactly one argument (an RHS-style expression listing",
            "covariates, e.g., `%s(x1 + x2 + x1:x2)`)."
          ),
          slot_name, fname, length(s) - 1L, fname
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_name, fname = fname,
                    nargs = length(s) - 1L)
      )
    }
    slot_key <- paste0(fname, "_formula")
    if (!is.null(out[[slot_key]])) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' calls `%s()` more than once. Only one occurrence of",
            "each AMM wrapper is admitted per slot; list multiple covariates",
            "inside a single call (e.g., `%s(x1 + x2)`)."
          ),
          slot_name, fname, fname
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_name, fname = fname)
      )
    }
    interior <- s[[2L]]
    sub_formula <- stats::as.formula(
      paste("~", deparse(interior, width.cutoff = 500L)),
      env = env
    )
    .gdpar_check_no_intercept_suppression(
      sub_formula,
      sprintf("%s -> %s()", slot_name, fname)
    )
    out[[slot_key]] <- sub_formula
  }
  out
}

#' Internal: split the RHS expression of an AMM formula into summands
#'
#' Recursively flattens a tree of binary \code{+} calls into a list of
#' top-level summands. Aborts on any non-\code{+} binary operator
#' (e.g., \code{-}, \code{*}, \code{/}, \code{:}) at the top level,
#' since the AMM RHS contract only admits additive composition of
#' the wrappers \code{a()}, \code{b()}, \code{W()}. The interior of
#' \code{a()} / \code{b()} is opaque to this splitter; arbitrary
#' formula syntax is allowed there (it is handed verbatim to
#' \code{\link[stats]{model.matrix}} downstream).
#'
#' @param expr An R expression representing the RHS of a formula.
#' @param slot_name Character scalar for error messages.
#' @return A list of R expressions, one per top-level summand.
#' @keywords internal
#' @noRd
.gdpar_split_amm_summands <- function(expr, slot_name) {
  if (is.call(expr) && length(expr) == 3L &&
      identical(expr[[1L]], as.symbol("+"))) {
    left  <- .gdpar_split_amm_summands(expr[[2L]], slot_name)
    right <- .gdpar_split_amm_summands(expr[[3L]], slot_name)
    return(c(left, right))
  }
  if (is.call(expr) && length(expr) == 3L) {
    op <- as.character(expr[[1L]])
    if (op %in% c("-", "*", "/", ":", "^", "|")) {
      gdpar_abort(
        sprintf(
          paste(
            "Slot '%s' uses the binary operator `%s` at the top level of its",
            "right-hand side. AMM canoniza solo composicion aditiva (`+`)",
            "de los wrappers `a()`, `b()`, `W()`; las interacciones y",
            "transformaciones van dentro del wrapper (e.g., `a(x1:x2)`)."
          ),
          slot_name, op
        ),
        class = "gdpar_input_error",
        data = list(slot = slot_name, operator = op)
      )
    }
  }
  list(expr)
}

#' Internal: convert a gdpar_formula_set to a named list of amm_spec
#'
#' For each slot of the supplied \code{gdpar_formula_set}, parses the
#' right-hand side via \code{.gdpar_parse_amm_formula()} and constructs
#' an \code{\link{amm_spec}} that combines the parsed components with
#' the external \code{W_basis} argument. The result is a named list
#' (names match the slot names) suitable for the K-individual dispatch
#' in \code{\link{gdpar}}.
#'
#' If any slot declares a \code{W()} term but no \code{W_basis} object
#' is supplied via the external \code{W} argument of \code{\link{gdpar}},
#' the function aborts with \code{gdpar_input_error} (decision P3,
#' sub-phase 8.3.3).
#'
#' For K = 1 (single slot) this function returns a one-element named
#' list; the caller is responsible for unwrapping when routing to the
#' legacy single-amm_spec path to preserve bit-exact backward
#' compatibility with Block 7.
#'
#' @param fs An object of class \code{gdpar_formula_set}.
#' @param W_basis_arg An object of class \code{W_basis} or \code{NULL}.
#' @return A named list of \code{amm_spec} objects.
#' @keywords internal
#' @noRd
.gdpar_formula_set_to_amm_spec_list <- function(fs, W_basis_arg = NULL) {
  if (!inherits(fs, "gdpar_formula_set")) {
    gdpar_abort(
      "Internal error: .gdpar_formula_set_to_amm_spec_list expected a gdpar_formula_set.",
      class = "gdpar_internal_error"
    )
  }
  out <- vector("list", length(fs$param_names))
  names(out) <- fs$param_names
  for (k in fs$param_names) {
    f <- fs$formulas[[k]]
    rhs_only <- if (length(f) == 3L) f[c(1L, 3L)] else f
    parsed <- .gdpar_parse_amm_formula(rhs_only, slot_name = k)
    W_for_slot <- NULL
    if (isTRUE(parsed$W_present)) {
      if (is.null(W_basis_arg)) {
        gdpar_abort(
          sprintf(
            paste(
              "Slot '%s' of the formula set declares `W()` but the external",
              "`W` argument of gdpar() is NULL. Pass an explicit",
              "`W = W_basis(...)` or remove the `W()` term from the slot."
            ),
            k
          ),
          class = "gdpar_input_error",
          data = list(slot = k)
        )
      }
      if (!inherits(W_basis_arg, "W_basis")) {
        gdpar_abort(
          paste(
            "The external `W` argument of gdpar() must be an object of class",
            "`W_basis` (see W_basis())."
          ),
          class = "gdpar_input_error",
          data = list(received = class(W_basis_arg))
        )
      }
      W_for_slot <- W_basis_arg
    }
    out[[k]] <- amm_spec(
      a = parsed$a_formula,
      b = parsed$b_formula,
      W = W_for_slot
    )
  }
  out
}

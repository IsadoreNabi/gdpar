#' Internal: assert that an object inherits a class
#'
#' @param x Object to check.
#' @param cls Character vector of acceptable classes.
#' @param arg_name Character scalar with the argument name as it
#'   appears in the calling function's signature, used in the error
#'   message.
#' @return Invisibly TRUE on success; otherwise signals a
#'   \code{gdpar_input_error}.
#' @keywords internal
#' @noRd
assert_inherits <- function(x, cls, arg_name) {
  if (!inherits(x, cls)) {
    gdpar_abort(
      sprintf(
        "Argument '%s' must inherit from class %s; received an object of class %s.",
        arg_name,
        paste(sQuote(cls), collapse = " or "),
        paste(sQuote(class(x)), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, expected = cls, received = class(x))
    )
  }
  invisible(TRUE)
}

#' Internal: assert that x is a single numeric value within a range
#'
#' @param x Object to check.
#' @param arg_name Character scalar with the argument name.
#' @param lower Numeric scalar with the inclusive lower bound, or
#'   \code{-Inf}.
#' @param upper Numeric scalar with the inclusive upper bound, or
#'   \code{Inf}.
#' @param allow_null Logical. If TRUE, \code{NULL} is accepted.
#' @return Invisibly TRUE on success.
#' @keywords internal
#' @noRd
assert_numeric_scalar <- function(x, arg_name, lower = -Inf, upper = Inf,
                                  allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) {
      return(invisible(TRUE))
    }
    gdpar_abort(
      sprintf("Argument '%s' must be a numeric scalar; received NULL.", arg_name),
      class = "gdpar_input_error",
      data = list(argument = arg_name)
    )
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    gdpar_abort(
      sprintf(
        "Argument '%s' must be a finite numeric scalar.", arg_name
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, received = x)
    )
  }
  if (x < lower || x > upper) {
    gdpar_abort(
      sprintf(
        "Argument '%s' must lie in [%g, %g]; received %g.",
        arg_name, lower, upper, x
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, lower = lower, upper = upper,
                  received = x)
    )
  }
  invisible(TRUE)
}

#' Internal: assert that x is a positive integer scalar
#'
#' @param x Object to check.
#' @param arg_name Character scalar with the argument name.
#' @param max Integer upper bound or \code{Inf}.
#' @return Invisibly TRUE on success.
#' @keywords internal
#' @noRd
assert_count <- function(x, arg_name, max = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != as.integer(x)) {
    gdpar_abort(
      sprintf(
        "Argument '%s' must be a positive integer scalar.", arg_name
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, received = x)
    )
  }
  if (x > max) {
    gdpar_abort(
      sprintf(
        "Argument '%s' must be at most %g; received %g.",
        arg_name, max, x
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, max = max, received = x)
    )
  }
  invisible(TRUE)
}

#' Internal: assert that x is a one-sided formula
#'
#' @param x Object to check.
#' @param arg_name Character scalar with the argument name.
#' @param allow_null Logical. If TRUE, \code{NULL} is accepted.
#' @return Invisibly TRUE on success.
#' @keywords internal
#' @noRd
assert_one_sided_formula <- function(x, arg_name, allow_null = TRUE) {
  if (is.null(x)) {
    if (allow_null) {
      return(invisible(TRUE))
    }
    gdpar_abort(
      sprintf("Argument '%s' must be a one-sided formula.", arg_name),
      class = "gdpar_input_error", data = list(argument = arg_name)
    )
  }
  if (!inherits(x, "formula") || length(x) != 2L) {
    gdpar_abort(
      sprintf(
        "Argument '%s' must be a one-sided formula such as ~ x1 + x2.",
        arg_name
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, received = deparse(x))
    )
  }
  invisible(TRUE)
}

#' Internal: assert that x is a data frame containing all variables
#'
#' @param x Object to check.
#' @param arg_name Character scalar with the argument name.
#' @param required_vars Character vector of variable names that must
#'   appear among the columns of x.
#' @return Invisibly TRUE on success.
#' @keywords internal
#' @noRd
assert_data_frame <- function(x, arg_name, required_vars = character()) {
  if (!is.data.frame(x)) {
    gdpar_abort(
      sprintf("Argument '%s' must be a data frame.", arg_name),
      class = "gdpar_input_error",
      data = list(argument = arg_name, received = class(x))
    )
  }
  missing_vars <- setdiff(required_vars, colnames(x))
  if (length(missing_vars) > 0) {
    gdpar_abort(
      sprintf(
        "Data frame '%s' is missing the following variable(s): %s.",
        arg_name, paste(sQuote(missing_vars), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(argument = arg_name, missing = missing_vars)
    )
  }
  invisible(TRUE)
}

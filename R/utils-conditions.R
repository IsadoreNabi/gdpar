#' Internal: signal a structured error
#'
#' Wraps \code{stop} with a class hierarchy so user code can catch
#' specific gdpar errors via \code{tryCatch}.
#'
#' @param message Character string with the error message.
#' @param class Character vector of class names appended before
#'   \code{c("gdpar_error", "error", "condition")}.
#' @param data Optional named list with additional fields embedded
#'   in the condition object for programmatic access.
#'
#' @return Never returns; signals a condition.
#' @keywords internal
#' @noRd
gdpar_abort <- function(message, class = character(), data = list()) {
  cond <- structure(
    class = c(class, "gdpar_error", "error", "condition"),
    list(message = message, call = sys.call(-1))
  )
  for (nm in names(data)) cond[[nm]] <- data[[nm]]
  stop(cond)
}

#' Internal: signal a structured warning
#'
#' Wraps \code{warning} with a class hierarchy so user code can catch
#' specific gdpar warnings via \code{withCallingHandlers}.
#'
#' @inheritParams gdpar_abort
#' @return Invisibly NULL; signals a condition.
#' @keywords internal
#' @noRd
gdpar_warn <- function(message, class = character(), data = list()) {
  cond <- structure(
    class = c(class, "gdpar_warning", "warning", "condition"),
    list(message = message, call = sys.call(-1))
  )
  for (nm in names(data)) cond[[nm]] <- data[[nm]]
  warning(cond)
}

#' Internal: signal a structured informative message
#'
#' Wraps \code{message} with a class hierarchy so user code can suppress
#' specific gdpar messages via \code{suppressMessages} or
#' \code{withCallingHandlers}.
#'
#' @inheritParams gdpar_abort
#' @return Invisibly NULL; signals a condition.
#' @keywords internal
#' @noRd
gdpar_inform <- function(message, class = character(), data = list()) {
  cond <- structure(
    class = c(class, "gdpar_message", "message", "condition"),
    list(message = paste0(message, "\n"), call = sys.call(-1))
  )
  for (nm in names(data)) cond[[nm]] <- data[[nm]]
  message(cond)
}

#' Internal: check that a Suggests-listed package is installed
#'
#' Used to enforce CRAN policy that Suggests packages may only be
#' invoked through conditional loading.
#'
#' @param pkg Character scalar with the package name.
#' @param reason Character scalar describing what the package is needed
#'   for. Embedded in the error message shown to the user.
#'
#' @return Invisibly TRUE if the package is available; otherwise
#'   signals a \code{gdpar_missing_dependency_error}.
#' @keywords internal
#' @noRd
require_suggested <- function(pkg, reason) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    gdpar_abort(
      sprintf(
        "Package '%s' is required to %s but is not installed. Please install it.",
        pkg, reason
      ),
      class = "gdpar_missing_dependency_error",
      data = list(package = pkg, reason = reason)
    )
  }
  invisible(TRUE)
}

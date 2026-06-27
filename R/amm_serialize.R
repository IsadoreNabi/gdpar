#' Serialize an amm_spec to a canonical plain-text file
#'
#' Write an \code{\link{amm_spec}} object to a text file in a canonical,
#' human-readable format suitable for version control, archival, and
#' bit-exact reproducibility. The file uses a small \code{key: value}
#' grammar with a mandatory version header and a hierarchical naming
#' convention for the per-dimension entries on the multivariate path.
#' The file is parsed by \code{\link{amm_load_spec}} via a dedicated
#' parser (no \code{source} or \code{eval} of the file contents), so
#' the serialized form is safe to load from untrusted locations.
#'
#' @param spec An object of class \code{amm_spec}.
#' @param path Character scalar giving the destination file path.
#'
#' @return Invisibly returns \code{path}.
#'
#' @details
#' The serialized format records exactly the constructor inputs of the
#' specification, not derived state. In particular, when the modulating
#' basis \code{W} has been materialised at a specific
#' \eqn{\theta_{\mathrm{ref}}} dimension (via the internal
#' \code{materialize_W_basis}), the materialised fields are not written;
#' the reconstructed object after \code{\link{amm_load_spec}} is the
#' unmaterialised \code{W_basis} corresponding to the same constructor
#' arguments, which is the form normally produced by
#' \code{\link{amm_spec}}.
#'
#' User-defined \code{W_basis} objects (\code{type = "user"}) cannot be
#' serialised because the evaluator is an arbitrary R function whose
#' definition cannot be canonised into the file format. An attempt to
#' serialise such a spec aborts with an informative error.
#'
#' @section File format:
#' The first non-empty line must be the version header
#' \code{# gdpar_spec_version: <version>}. Subsequent lines are either
#' comments starting with \code{#}, empty lines (ignored), or records of
#' the form \code{key: value}. Recognised keys and value grammars:
#' \itemize{
#'   \item \code{p}: positive integer.
#'   \item \code{a}, \code{b}: \code{NULL} or a one-sided formula
#'     literal such as \code{~ x1 + x2}. Used on the scalar path; set
#'     to \code{NULL} on the multivariate path.
#'   \item \code{x_vars}: \code{NULL} or a literal of the form
#'     \code{c("x1", "x2", ...)}.
#'   \item \code{W.type}: \code{NULL}, \code{polynomial} or
#'     \code{bspline}. \code{user} is not serialisable.
#'   \item \code{W.degree}: positive integer (required for
#'     \code{polynomial} and \code{bspline}).
#'   \item \code{W.knots}: \code{c(...)} of numerics (bspline with
#'     interior knots only).
#'   \item \code{W.df}: positive integer (bspline with df only).
#'   \item \code{dims.K.a}, \code{dims.K.b}: same grammar as \code{a},
#'     \code{b}, for \code{K} in \code{1:p}. Required when \code{p > 1}.
#' }
#'
#' @section Version policy:
#' The version field is checked strictly against the running package
#' version. Until the package reaches its first stable release a
#' mismatch is treated as an error; a subsequent release will introduce
#' a forward-compatible upgrade path.
#'
#' @examples
#' \donttest{
#' spec <- amm_spec(a = ~ x1 + x2, b = ~ x1,
#'                  W = W_basis(type = "polynomial", degree = 2))
#' tmp  <- tempfile(fileext = ".gdpar")
#' amm_save_spec(spec, tmp)
#' spec2 <- amm_load_spec(tmp)
#' identical(spec$level, spec2$level)
#' }
#'
#' @seealso \code{\link{amm_load_spec}}, \code{\link{amm_spec}},
#'   \code{\link{W_basis}}
#' @export
amm_save_spec <- function(spec, path) {
  assert_inherits(spec, "amm_spec", "spec")
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    gdpar_abort(
      "Argument 'path' must be a non-empty character scalar.",
      class = "gdpar_input_error",
      data = list(argument = "path", received = path)
    )
  }
  W <- spec[["W"]]
  if (!is.null(W) && identical(W[["type"]], "user")) {
    gdpar_abort(
      paste(
        "User-defined W_basis (type = 'user') cannot be serialised:",
        "the evaluator is an arbitrary R function. Use 'polynomial' or",
        "'bspline' bases, or reconstruct the user evaluator manually",
        "after loading a scaffold spec."
      ),
      class = "gdpar_input_error"
    )
  }
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  lines <- c(paste0("# gdpar_spec_version: ", pkg_ver),
             paste0("p: ", spec[["p"]]))
  if (isTRUE(spec[["p"]] == 1L)) {
    lines <- c(lines,
               paste0("a: ", .serialize_formula(spec[["a"]])),
               paste0("b: ", .serialize_formula(spec[["b"]])))
  } else {
    lines <- c(lines, "a: NULL", "b: NULL")
  }
  lines <- c(lines,
             paste0("x_vars: ", .serialize_char_vec(spec[["x_vars"]])))
  if (is.null(W)) {
    lines <- c(lines, "W.type: NULL")
  } else {
    lines <- c(lines, paste0("W.type: ", W[["type"]]))
    if (identical(W[["type"]], "polynomial")) {
      lines <- c(lines, paste0("W.degree: ", as.integer(W[["degree"]])))
    } else if (identical(W[["type"]], "bspline")) {
      lines <- c(lines, paste0("W.degree: ", as.integer(W[["degree"]])))
      if (!is.null(W[["knots"]])) {
        lines <- c(lines,
                   paste0("W.knots: ", .serialize_num_vec(W[["knots"]])))
      }
      if (!is.null(W[["df"]])) {
        lines <- c(lines, paste0("W.df: ", as.integer(W[["df"]])))
      }
    }
  }
  if (isTRUE(spec[["p"]] > 1L)) {
    for (k in seq_len(spec[["p"]])) {
      entry <- spec[["dims"]][[k]]
      lines <- c(lines,
                 paste0("dims.", k, ".a: ",
                        .serialize_formula(entry[["a"]])),
                 paste0("dims.", k, ".b: ",
                        .serialize_formula(entry[["b"]])))
    }
  }
  writeLines(lines, con = path)
  invisible(path)
}

#' Load a canonical amm_spec file produced by amm_save_spec
#'
#' Parse a file written by \code{\link{amm_save_spec}} and reconstruct
#' the corresponding \code{\link{amm_spec}} object. Parsing is purely
#' lexical (no \code{source} or \code{eval} of the file contents) and
#' validates the version header strictly against the running package
#' version.
#'
#' @param path Character scalar giving the file path.
#'
#' @return An object of class \code{amm_spec}.
#'
#' @section Errors:
#' Aborts with a \code{gdpar_input_error} for: missing or unreadable
#' file; absent or malformed version header; version mismatch; malformed
#' record line (missing colon); unknown or missing required keys;
#' invalid value grammar for a recognised key; \code{W.type = "user"}
#' (not serialisable in the canonical form); \code{dims.K.*} records
#' outside \code{1:p}.
#'
#' @examples
#' \donttest{
#' spec <- amm_spec(a = ~ x1 + x2,
#'                  W = W_basis(type = "polynomial", degree = 2))
#' tmp  <- tempfile(fileext = ".gdpar")
#' amm_save_spec(spec, tmp)
#' spec2 <- amm_load_spec(tmp)
#' print(spec2)
#' }
#'
#' @seealso \code{\link{amm_save_spec}}, \code{\link{amm_spec}}
#' @export
amm_load_spec <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    gdpar_abort(
      "Argument 'path' must be a non-empty character scalar.",
      class = "gdpar_input_error",
      data = list(argument = "path", received = path)
    )
  }
  if (!file.exists(path)) {
    gdpar_abort(
      sprintf("File not found: %s.", path),
      class = "gdpar_input_error",
      data = list(path = path)
    )
  }
  raw_lines <- readLines(path, warn = FALSE)
  version_idx <- which(grepl("^\\s*#\\s*gdpar_spec_version\\s*:", raw_lines))
  if (length(version_idx) == 0L) {
    gdpar_abort(
      paste(
        "Missing required header '# gdpar_spec_version: <version>'.",
        "The file does not appear to be a canonical gdpar spec."
      ),
      class = "gdpar_input_error"
    )
  }
  ver_value <- trimws(sub("^\\s*#\\s*gdpar_spec_version\\s*:",
                          "", raw_lines[[version_idx[[1L]]]]))
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  if (!identical(ver_value, pkg_ver)) {
    gdpar_abort(
      sprintf(paste(
        "Version mismatch: file declares gdpar %s but the loaded package",
        "is %s. Bit-exact reproducibility across versions is not",
        "guaranteed for this development release."),
        ver_value, pkg_ver),
      class = "gdpar_input_error",
      data = list(file_version = ver_value, package_version = pkg_ver)
    )
  }
  records <- list()
  for (i in seq_along(raw_lines)) {
    raw <- raw_lines[[i]]
    s <- trimws(raw)
    if (!nzchar(s)) next
    if (startsWith(s, "#")) next
    pos <- regexpr(":", raw, fixed = TRUE)
    if (pos == -1L) {
      gdpar_abort(
        sprintf("Malformed record at line %d (no ':'): %s", i, raw),
        class = "gdpar_input_error",
        data = list(line = i, raw = raw)
      )
    }
    key <- trimws(substr(raw, 1L, pos - 1L))
    val <- trimws(substr(raw, pos + 1L, nchar(raw)))
    if (!nzchar(key)) {
      gdpar_abort(
        sprintf("Empty key at line %d.", i),
        class = "gdpar_input_error"
      )
    }
    if (!is.null(records[[key]])) {
      gdpar_abort(
        sprintf("Duplicate key '%s' at line %d (previously seen at line %d).",
                key, i, records[[key]][["line"]]),
        class = "gdpar_input_error"
      )
    }
    records[[key]] <- list(value = val, line = i)
  }
  recognised_prefixes <- c("p", "a", "b", "x_vars",
                           "W.type", "W.degree", "W.knots", "W.df")
  for (key in names(records)) {
    if (key %in% recognised_prefixes) next
    if (grepl("^dims\\.[0-9]+\\.[ab]$", key)) next
    gdpar_abort(
      sprintf("Unknown key '%s' at line %d.",
              key, records[[key]][["line"]]),
      class = "gdpar_input_error",
      data = list(key = key)
    )
  }
  if (is.null(records[["p"]])) {
    gdpar_abort(
      "Required key 'p' not found.",
      class = "gdpar_input_error"
    )
  }
  p_val <- .parse_int(records[["p"]][["value"]], "p",
                      records[["p"]][["line"]])
  if (p_val < 1L) {
    gdpar_abort(
      sprintf("Key 'p' must be >= 1; got %d at line %d.",
              p_val, records[["p"]][["line"]]),
      class = "gdpar_input_error"
    )
  }
  a_rec <- records[["a"]]
  b_rec <- records[["b"]]
  if (is.null(a_rec)) {
    gdpar_abort("Required key 'a' not found.", class = "gdpar_input_error")
  }
  if (is.null(b_rec)) {
    gdpar_abort("Required key 'b' not found.", class = "gdpar_input_error")
  }
  a_scalar <- .parse_formula(a_rec[["value"]], "a", a_rec[["line"]])
  b_scalar <- .parse_formula(b_rec[["value"]], "b", b_rec[["line"]])
  xv_rec <- records[["x_vars"]]
  x_vars <- if (is.null(xv_rec)) NULL else {
    .parse_char_vec(xv_rec[["value"]], "x_vars", xv_rec[["line"]])
  }
  W <- .parse_W_records(records)
  if (p_val == 1L) {
    dim_keys <- grep("^dims\\.", names(records), value = TRUE)
    if (length(dim_keys) > 0L) {
      gdpar_abort(
        sprintf(paste(
          "Key 'p' is 1 but the file contains multivariate dims records:",
          "%s."),
          paste(sQuote(dim_keys), collapse = ", ")),
        class = "gdpar_input_error"
      )
    }
    return(amm_spec(a = a_scalar, b = b_scalar, W = W,
                    x_vars = x_vars, p = 1L))
  }
  if (!is.null(a_scalar) || !is.null(b_scalar)) {
    gdpar_abort(
      paste(
        "Keys 'a' and 'b' must be NULL when p > 1; per-dimension values",
        "go in dims.K.a and dims.K.b records."
      ),
      class = "gdpar_input_error"
    )
  }
  dims_list <- vector("list", p_val)
  for (k in seq_len(p_val)) {
    a_key <- sprintf("dims.%d.a", k)
    b_key <- sprintf("dims.%d.b", k)
    a_rec_k <- records[[a_key]]
    b_rec_k <- records[[b_key]]
    if (is.null(a_rec_k)) {
      gdpar_abort(
        sprintf("Required key '%s' not found for p = %d.", a_key, p_val),
        class = "gdpar_input_error"
      )
    }
    if (is.null(b_rec_k)) {
      gdpar_abort(
        sprintf("Required key '%s' not found for p = %d.", b_key, p_val),
        class = "gdpar_input_error"
      )
    }
    a_k <- .parse_formula(a_rec_k[["value"]], a_key, a_rec_k[["line"]])
    b_k <- .parse_formula(b_rec_k[["value"]], b_key, b_rec_k[["line"]])
    dims_list[[k]] <- list(a = a_k, b = b_k)
  }
  extra_dims <- grep("^dims\\.", names(records), value = TRUE)
  extra_dims <- extra_dims[!grepl(sprintf("^dims\\.[1-%d]\\.[ab]$|^dims\\.%d\\.[ab]$",
                                          p_val - 1L, p_val),
                                  extra_dims)]
  if (length(extra_dims) > 0L) {
    bad <- character(0)
    for (key in extra_dims) {
      m <- regmatches(key, regexec("^dims\\.([0-9]+)\\.[ab]$", key))[[1L]]
      if (length(m) >= 2L) {
        k_idx <- as.integer(m[[2L]])
        if (is.na(k_idx) || k_idx < 1L || k_idx > p_val) {
          bad <- c(bad, key)
        }
      } else {
        bad <- c(bad, key)
      }
    }
    if (length(bad) > 0L) {
      gdpar_abort(
        sprintf("dims.K.* records outside 1:%d: %s.",
                p_val, paste(sQuote(bad), collapse = ", ")),
        class = "gdpar_input_error"
      )
    }
  }
  amm_spec(W = W, x_vars = x_vars, p = p_val, dims = dims_list)
}

.serialize_formula <- function(f) {
  if (is.null(f)) return("NULL")
  paste(deparse(f, width.cutoff = 500L), collapse = " ")
}

.parse_formula <- function(value, key, line_no) {
  value <- trimws(value)
  if (identical(value, "NULL")) return(NULL)
  if (!startsWith(value, "~")) {
    gdpar_abort(
      sprintf(paste(
        "Value for key '%s' must be NULL or a one-sided formula starting",
        "with '~'; got %s at line %d."),
        key, sQuote(value), line_no),
      class = "gdpar_input_error",
      data = list(key = key, value = value, line = line_no)
    )
  }
  out <- tryCatch(
    stats::as.formula(value),
    error = function(e) {
      gdpar_abort(
        sprintf(paste(
          "Could not parse the value of key '%s' at line %d as a formula:",
          "%s"),
          key, line_no, conditionMessage(e)),
        class = "gdpar_input_error"
      )
    }
  )
  if (length(out) != 2L) {
    gdpar_abort(
      sprintf(paste(
        "Value for key '%s' must be a one-sided formula; got a two-sided",
        "expression at line %d."),
        key, line_no),
      class = "gdpar_input_error"
    )
  }
  out
}

.serialize_char_vec <- function(x) {
  if (is.null(x)) return("NULL")
  if (any(grepl("[\"\\\\\n\r]", x))) {
    gdpar_abort(
      paste(
        "Cannot serialise x_vars containing double quotes, backslashes,",
        "or newlines."
      ),
      class = "gdpar_input_error"
    )
  }
  paste0("c(", paste(sprintf("\"%s\"", x), collapse = ", "), ")")
}

.parse_char_vec <- function(value, key, line_no) {
  value <- trimws(value)
  if (identical(value, "NULL")) return(NULL)
  if (!grepl("^c\\(.*\\)$", value)) {
    gdpar_abort(
      sprintf(paste(
        "Value for key '%s' must be NULL or a c(...) literal at line",
        "%d; got %s."),
        key, line_no, sQuote(value)),
      class = "gdpar_input_error"
    )
  }
  inner <- sub("^c\\((.*)\\)$", "\\1", value)
  matches <- regmatches(inner, gregexpr("\"([^\"]*)\"", inner))[[1L]]
  if (length(matches) == 0L) {
    if (!nzchar(trimws(inner))) return(character(0))
    gdpar_abort(
      sprintf(paste(
        "Could not parse string entries of '%s' at line %d:",
        "expected quoted tokens."),
        key, line_no),
      class = "gdpar_input_error"
    )
  }
  gsub("\"", "", matches, fixed = TRUE)
}

.serialize_num_vec <- function(x) {
  if (is.null(x)) return("NULL")
  if (length(x) == 0L) return("c()")
  formatted <- vapply(x, function(v) sprintf("%.17g", v), character(1L))
  paste0("c(", paste(formatted, collapse = ", "), ")")
}

.parse_num_vec <- function(value, key, line_no) {
  value <- trimws(value)
  if (identical(value, "NULL")) return(NULL)
  if (!grepl("^c\\(.*\\)$", value)) {
    gdpar_abort(
      sprintf(paste(
        "Value for key '%s' must be NULL or a c(...) literal at line",
        "%d; got %s."),
        key, line_no, sQuote(value)),
      class = "gdpar_input_error"
    )
  }
  inner <- trimws(sub("^c\\((.*)\\)$", "\\1", value))
  if (!nzchar(inner)) return(numeric(0))
  parts <- trimws(strsplit(inner, ",", fixed = TRUE)[[1L]])
  vals <- suppressWarnings(as.numeric(parts))
  if (any(is.na(vals))) {
    gdpar_abort(
      sprintf(paste(
        "Could not parse numeric entries of '%s' at line %d:",
        "non-numeric tokens present."),
        key, line_no),
      class = "gdpar_input_error"
    )
  }
  vals
}

.parse_int <- function(value, key, line_no) {
  value <- trimws(value)
  v <- suppressWarnings(as.numeric(value))
  if (is.na(v) || !is.finite(v) || v != as.integer(v)) {
    gdpar_abort(
      sprintf("Value for key '%s' must be an integer at line %d; got %s.",
              key, line_no, sQuote(value)),
      class = "gdpar_input_error"
    )
  }
  as.integer(v)
}

.parse_W_records <- function(records) {
  Wt_rec <- records[["W.type"]]
  if (is.null(Wt_rec) || identical(Wt_rec[["value"]], "NULL")) {
    return(NULL)
  }
  W_type <- Wt_rec[["value"]]
  if (identical(W_type, "user")) {
    gdpar_abort(
      paste(
        "W.type = 'user' is not supported in the canonical file format:",
        "user-defined bases reference arbitrary R functions."
      ),
      class = "gdpar_input_error"
    )
  }
  if (!W_type %in% c("polynomial", "bspline")) {
    gdpar_abort(
      sprintf(paste(
        "Unsupported W.type %s at line %d. Supported types in the",
        "canonical format: polynomial, bspline."),
        sQuote(W_type), Wt_rec[["line"]]),
      class = "gdpar_input_error"
    )
  }
  deg_rec <- records[["W.degree"]]
  if (is.null(deg_rec)) {
    gdpar_abort(
      sprintf("Key 'W.degree' is required when W.type = '%s'.", W_type),
      class = "gdpar_input_error"
    )
  }
  W_degree <- .parse_int(deg_rec[["value"]], "W.degree", deg_rec[["line"]])
  if (identical(W_type, "polynomial")) {
    return(W_basis(type = "polynomial", degree = W_degree))
  }
  knots_rec <- records[["W.knots"]]
  df_rec <- records[["W.df"]]
  if (!is.null(knots_rec) && !is.null(df_rec)) {
    gdpar_abort(
      paste(
        "Keys 'W.knots' and 'W.df' are mutually exclusive for",
        "W.type = 'bspline'."
      ),
      class = "gdpar_input_error"
    )
  }
  if (!is.null(knots_rec)) {
    knots <- .parse_num_vec(knots_rec[["value"]], "W.knots",
                            knots_rec[["line"]])
    return(W_basis(type = "bspline", degree = W_degree, knots = knots))
  }
  if (!is.null(df_rec)) {
    df_val <- .parse_int(df_rec[["value"]], "W.df", df_rec[["line"]])
    return(W_basis(type = "bspline", degree = W_degree, df = df_val))
  }
  gdpar_abort(
    paste(
      "For W.type = 'bspline' one of 'W.knots' or 'W.df' must be",
      "supplied."
    ),
    class = "gdpar_input_error"
  )
}

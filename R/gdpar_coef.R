#' Internal: build a per-component data.frame of posterior summaries
#'
#' Given a matrix of draws (rows = samples, cols = parameters) and a
#' character vector of column labels, return a tidy data.frame with one
#' row per column and the standard summary columns
#' (\code{mean}, \code{q05}, \code{q50}, \code{q95}). Used by the multi
#' and scalar branches of \code{coef.gdpar_fit()}.
#'
#' @keywords internal
#' @noRd
build_coef_term_df <- function(draws_mat, term_names) {
  if (is.null(draws_mat) || ncol(draws_mat) == 0L) {
    return(NULL)
  }
  if (length(term_names) != ncol(draws_mat)) {
    gdpar_abort(
      sprintf(
        "Internal error: term_names length (%d) does not match draws columns (%d).",
        length(term_names), ncol(draws_mat)
      ),
      class = "gdpar_internal_error"
    )
  }
  means <- colMeans(draws_mat)
  qs <- apply(draws_mat, 2L, stats::quantile,
              probs = c(0.05, 0.50, 0.95), names = FALSE,
              na.rm = TRUE)
  data.frame(
    term  = term_names,
    mean  = unname(means),
    q05   = unname(qs[1L, ]),
    q50   = unname(qs[2L, ]),
    q95   = unname(qs[3L, ]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Internal: build a tidy data.frame for a W block (per coordinate)
#'
#' For a draws matrix of dimensions \code{(S, dim_W_per_k * d_x)} that
#' represents the row-major flattening of a \code{W_per_k_dim x d_x}
#' matrix, return a long-tidy data.frame with cols
#' \code{(basis_idx, x_name, mean, q05, q50, q95)}.
#'
#' @keywords internal
#' @noRd
build_coef_W_df <- function(draws_mat, basis_dim, x_names) {
  if (is.null(draws_mat) || ncol(draws_mat) == 0L ||
      basis_dim == 0L || length(x_names) == 0L) {
    return(NULL)
  }
  if (ncol(draws_mat) != basis_dim * length(x_names)) {
    gdpar_abort(
      sprintf(
        "Internal error: W draws columns (%d) != basis_dim * d_x (%d * %d).",
        ncol(draws_mat), basis_dim, length(x_names)
      ),
      class = "gdpar_internal_error"
    )
  }
  means <- colMeans(draws_mat)
  qs <- apply(draws_mat, 2L, stats::quantile,
              probs = c(0.05, 0.50, 0.95), names = FALSE,
              na.rm = TRUE)
  basis_idx <- rep(seq_len(basis_dim), times = length(x_names))
  x_name <- rep(x_names, each = basis_dim)
  data.frame(
    basis_idx = basis_idx,
    x_name    = x_name,
    mean      = unname(means),
    q05       = unname(qs[1L, ]),
    q50       = unname(qs[2L, ]),
    q95       = unname(qs[3L, ]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Internal: build the theta_ref summary data.frame
#'
#' For a draws matrix of dimensions \code{(S, p)} return a tidy
#' data.frame with one row per coordinate and cols
#' \code{(k, mean, q05, q50, q95)}.
#'
#' @keywords internal
#' @noRd
build_coef_theta_ref_df <- function(draws_mat, p) {
  if (ncol(draws_mat) != p) {
    gdpar_abort(
      sprintf(
        "Internal error: theta_ref draws columns (%d) != p (%d).",
        ncol(draws_mat), p
      ),
      class = "gdpar_internal_error"
    )
  }
  means <- colMeans(draws_mat)
  qs <- apply(draws_mat, 2L, stats::quantile,
              probs = c(0.05, 0.50, 0.95), names = FALSE,
              na.rm = TRUE)
  data.frame(
    k     = seq_len(p),
    mean  = unname(means),
    q05   = unname(qs[1L, ]),
    q50   = unname(qs[2L, ]),
    q95   = unname(qs[3L, ]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Internal: build the theta_ref summary data.frame under grouping
#'
#' Like build_coef_theta_ref_df but accepts the per-group axis of Block
#' 6.5. The input is either an (S, J_groups) matrix (univariate p = 1)
#' or an (S, J_groups, p) array (multivariate). The output is a tidy
#' data.frame with cols (g, k, mean, q05, q50, q95) when J_groups > 1
#' and cols (k, mean, q05, q50, q95) when J_groups == 1 (matching the
#' Block 6 backward-compat schema bit-exactly).
#'
#' @keywords internal
#' @noRd
build_coef_theta_ref_df_grouped <- function(arr, J_groups, p) {
  if (J_groups == 1L) {
    if (is.array(arr) && length(dim(arr)) == 3L) {
      flat <- matrix(arr[, 1L, ], nrow = dim(arr)[1L], ncol = p)
    } else if (is.matrix(arr) && ncol(arr) == 1L && p == 1L) {
      flat <- arr
    } else if (is.matrix(arr) && ncol(arr) == p) {
      flat <- arr
    } else {
      gdpar_abort(
        "Internal error: unexpected shape for theta_ref draws under J_groups = 1.",
        class = "gdpar_internal_error"
      )
    }
    return(build_coef_theta_ref_df(flat, p))
  }
  if (is.matrix(arr)) {
    if (ncol(arr) != J_groups || p != 1L) {
      gdpar_abort(
        "Internal error: theta_ref matrix shape does not match (J_groups, 1).",
        class = "gdpar_internal_error"
      )
    }
    means <- colMeans(arr)
    qs <- apply(arr, 2L, stats::quantile,
                probs = c(0.05, 0.50, 0.95), names = FALSE,
                na.rm = TRUE)
    return(data.frame(
      g     = seq_len(J_groups),
      k     = rep(1L, J_groups),
      mean  = unname(means),
      q05   = unname(qs[1L, ]),
      q50   = unname(qs[2L, ]),
      q95   = unname(qs[3L, ]),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }
  if (!is.array(arr) || length(dim(arr)) != 3L ||
      dim(arr)[2L] != J_groups || dim(arr)[3L] != p) {
    gdpar_abort(
      sprintf(
        "Internal error: theta_ref array dims do not match (S, %d, %d).",
        J_groups, p
      ),
      class = "gdpar_internal_error"
    )
  }
  out_rows <- list()
  for (g in seq_len(J_groups)) {
    for (k in seq_len(p)) {
      v <- arr[, g, k]
      qs <- stats::quantile(v, probs = c(0.05, 0.50, 0.95),
                            names = FALSE, na.rm = TRUE)
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        g     = g,
        k     = k,
        mean  = mean(v, na.rm = TRUE),
        q05   = qs[1L],
        q50   = qs[2L],
        q95   = qs[3L],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  do.call(rbind, out_rows)
}

#' Internal: build a hyperparameter summary data.frame
#'
#' Common shape for mu_theta_ref and sigma_theta_ref: one row per coord
#' k, cols (k, mean, q05, q50, q95). Input is the (S, p) matrix returned
#' by .extract_mu_sigma_theta_ref.
#'
#' @keywords internal
#' @noRd
build_coef_hyper_df <- function(draws_mat, p) {
  if (is.null(draws_mat)) return(NULL)
  if (ncol(draws_mat) != p) {
    gdpar_abort(
      sprintf(
        "Internal error: hyper draws columns (%d) != p (%d).",
        ncol(draws_mat), p
      ),
      class = "gdpar_internal_error"
    )
  }
  means <- colMeans(draws_mat)
  qs <- apply(draws_mat, 2L, stats::quantile,
              probs = c(0.05, 0.50, 0.95), names = FALSE,
              na.rm = TRUE)
  data.frame(
    k     = seq_len(p),
    mean  = unname(means),
    q05   = unname(qs[1L, ]),
    q50   = unname(qs[2L, ]),
    q95   = unname(qs[3L, ]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Internal: validate a per-k coef list slot
#'
#' Checks that the slot is either \code{NULL} (component absent at the
#' AMM level) or a list of length \code{p} whose entries are
#' \code{NULL} (coord-inactive) or data.frames with the expected
#' columns.
#'
#' @keywords internal
#' @noRd
validate_coef_slot <- function(slot, p, expected_cols, slot_name) {
  if (is.null(slot)) {
    return(invisible(NULL))
  }
  if (!is.list(slot) || length(slot) != p) {
    gdpar_abort(
      sprintf(
        "Slot '%s' must be NULL or a list of length p (%d); got %s of length %d.",
        slot_name, p, class(slot)[1L], length(slot)
      ),
      class = "gdpar_internal_error"
    )
  }
  for (k in seq_len(p)) {
    elem <- slot[[k]]
    if (is.null(elem)) next
    if (!is.data.frame(elem)) {
      gdpar_abort(
        sprintf(
          "Slot '%s'[[%d]] must be NULL or data.frame; got %s.",
          slot_name, k, class(elem)[1L]
        ),
        class = "gdpar_internal_error"
      )
    }
    missing_cols <- setdiff(expected_cols, names(elem))
    if (length(missing_cols) > 0L) {
      gdpar_abort(
        sprintf(
          "Slot '%s'[[%d]] missing required cols: %s.",
          slot_name, k, paste(missing_cols, collapse = ", ")
        ),
        class = "gdpar_internal_error"
      )
    }
  }
  invisible(NULL)
}

#' Internal: validate the theta_ref data.frame of a gdpar_coef
#'
#' @keywords internal
#' @noRd
validate_coef_theta_ref <- function(theta_ref_df, p, J_groups = 1L) {
  if (!is.data.frame(theta_ref_df)) {
    gdpar_abort(
      "Slot 'theta_ref' must be a data.frame.",
      class = "gdpar_internal_error"
    )
  }
  has_group <- "g" %in% names(theta_ref_df)
  expected <- c("k", "mean", "q05", "q50", "q95")
  if (has_group) expected <- c("g", expected)
  missing_cols <- setdiff(expected, names(theta_ref_df))
  if (length(missing_cols) > 0L) {
    gdpar_abort(
      sprintf(
        "Slot 'theta_ref' missing required cols: %s.",
        paste(missing_cols, collapse = ", ")
      ),
      class = "gdpar_internal_error"
    )
  }
  expected_rows <- if (has_group) J_groups * p else p
  if (nrow(theta_ref_df) != expected_rows) {
    gdpar_abort(
      sprintf(
        "Slot 'theta_ref' has %d rows; expected %d (= J_groups * p with J_groups = %d, p = %d).",
        nrow(theta_ref_df), expected_rows, J_groups, p
      ),
      class = "gdpar_internal_error"
    )
  }
  if (!has_group) {
    if (!identical(sort(as.integer(theta_ref_df$k)), seq_len(p))) {
      gdpar_abort(
        "Slot 'theta_ref$k' must contain integers 1..p without gaps.",
        class = "gdpar_internal_error"
      )
    }
  } else {
    if (!identical(sort(unique(as.integer(theta_ref_df$g))),
                   seq_len(J_groups))) {
      gdpar_abort(
        "Slot 'theta_ref$g' must contain integers 1..J_groups without gaps.",
        class = "gdpar_internal_error"
      )
    }
    if (!identical(sort(unique(as.integer(theta_ref_df$k))),
                   seq_len(p))) {
      gdpar_abort(
        "Slot 'theta_ref$k' must contain integers 1..p without gaps.",
        class = "gdpar_internal_error"
      )
    }
  }
  invisible(NULL)
}

#' Internal: validate a hyperparameter data.frame (mu_theta_ref or sigma_theta_ref)
#'
#' Expected cols (k, mean, q05, q50, q95), p rows. NULL slots are
#' accepted (component inactive: no grouping).
#'
#' @keywords internal
#' @noRd
validate_coef_hyper <- function(df, p, slot_name) {
  if (is.null(df)) return(invisible(NULL))
  if (!is.data.frame(df)) {
    gdpar_abort(
      sprintf("Slot '%s' must be NULL or a data.frame.", slot_name),
      class = "gdpar_internal_error"
    )
  }
  expected <- c("k", "mean", "q05", "q50", "q95")
  missing_cols <- setdiff(expected, names(df))
  if (length(missing_cols) > 0L) {
    gdpar_abort(
      sprintf(
        "Slot '%s' missing required cols: %s.",
        slot_name, paste(missing_cols, collapse = ", ")
      ),
      class = "gdpar_internal_error"
    )
  }
  if (nrow(df) != p) {
    gdpar_abort(
      sprintf(
        "Slot '%s' has %d rows; expected p (%d).",
        slot_name, nrow(df), p
      ),
      class = "gdpar_internal_error"
    )
  }
  invisible(NULL)
}

#' Internal constructor for the gdpar_coef class
#'
#' Builds an object of class \code{gdpar_coef} from already-summarized
#' per-component data.frames. All validation is performed here so that
#' downstream methods can assume a well-formed object.
#'
#' @param theta_ref data.frame with cols
#'   \code{(k, mean, q05, q50, q95)}, p rows.
#' @param a NULL (component absent in the AMM) or a list of length p
#'   whose entries are NULL (coord inactive) or data.frames with cols
#'   \code{(term, mean, q05, q50, q95)}.
#' @param b Same conventions as \code{a}.
#' @param W NULL or a list of length p whose entries are NULL or
#'   data.frames with cols
#'   \code{(basis_idx, x_name, mean, q05, q50, q95)}.
#' @param p Integer scalar >= 1.
#' @param summary_stats Character vector of statistics included.
#' @return Object of class \code{gdpar_coef}.
#' @keywords internal
#' @noRd
new_gdpar_coef <- function(theta_ref, a = NULL, b = NULL, W = NULL,
                           p,
                           mu_theta_ref = NULL,
                           sigma_theta_ref = NULL,
                           J_groups = 1L,
                           group_levels = NULL,
                           summary_stats = c("mean", "q05", "q50", "q95")) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 1L ||
      p != as.integer(p)) {
    gdpar_abort(
      sprintf("'p' must be a positive integer scalar; got %s.",
              paste(deparse(p), collapse = " ")),
      class = "gdpar_internal_error"
    )
  }
  p <- as.integer(p)
  J_groups <- as.integer(J_groups)
  validate_coef_theta_ref(theta_ref, p, J_groups)
  validate_coef_slot(a, p, c("term", "mean", "q05", "q50", "q95"), "a")
  validate_coef_slot(b, p, c("term", "mean", "q05", "q50", "q95"), "b")
  validate_coef_slot(W, p,
                     c("basis_idx", "x_name", "mean", "q05", "q50", "q95"),
                     "W")
  validate_coef_hyper(mu_theta_ref, p, "mu_theta_ref")
  validate_coef_hyper(sigma_theta_ref, p, "sigma_theta_ref")
  structure(
    list(
      theta_ref       = theta_ref,
      a               = a,
      b               = b,
      W               = W,
      p               = p,
      mu_theta_ref    = mu_theta_ref,
      sigma_theta_ref = sigma_theta_ref,
      J_groups        = J_groups,
      group_levels    = group_levels,
      summary_stats   = summary_stats
    ),
    class = "gdpar_coef"
  )
}

#' Internal: count coordinates with an active component slot
#'
#' @keywords internal
#' @noRd
count_active_coords <- function(slot) {
  if (is.null(slot)) return(0L)
  sum(vapply(slot, function(e) {
    !is.null(e) && is.data.frame(e) && nrow(e) > 0L
  }, logical(1L)))
}

#' Print method for gdpar_coef objects
#'
#' Three verbosity levels: \code{"global"} (default) shows the
#' theta_ref summary plus active-component counts;
#' \code{"coord"} appends per-coordinate component means;
#' \code{"full"} appends every per-coordinate data.frame with all
#' summary statistics.
#'
#' @param x Object of class \code{gdpar_coef}.
#' @param level One of \code{"global"}, \code{"coord"}, \code{"full"}.
#' @param digits Integer scalar passed to \code{format()} for numeric
#'   columns.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_coef <- function(x,
                             level = c("global", "coord", "full"),
                             digits = 4L,
                             ...) {
  level <- match.arg(level)
  cat("<gdpar_coef>\n")
  cat("  p                 : ", x$p, "\n", sep = "")
  J_groups <- if (is.null(x$J_groups)) 1L else x$J_groups
  if (J_groups > 1L) {
    cat("  J_groups          : ", J_groups, "\n", sep = "")
  }
  cat("  summary_stats     : ",
      paste(x$summary_stats, collapse = ", "), "\n", sep = "")
  n_a <- count_active_coords(x$a)
  n_b <- count_active_coords(x$b)
  n_W <- count_active_coords(x$W)
  cat(sprintf(
    "  components active : a(%d/%d) b(%d/%d) W(%d/%d)\n",
    n_a, x$p, n_b, x$p, n_W, x$p
  ))
  if (!is.null(x$mu_theta_ref)) {
    cat("\n  mu_theta_ref:\n")
    print(format_coef_df(x$mu_theta_ref, digits = digits),
          row.names = FALSE)
  }
  if (!is.null(x$sigma_theta_ref)) {
    cat("\n  sigma_theta_ref:\n")
    print(format_coef_df(x$sigma_theta_ref, digits = digits),
          row.names = FALSE)
  }
  cat("\n  theta_ref:\n")
  print(format_coef_df(x$theta_ref, digits = digits), row.names = FALSE)
  if (level == "global") {
    cat("\n  use print(x, level = \"coord\") or print(x, level = \"full\")",
        " for per-coordinate detail; as.data.frame(x) for the long-tidy table.\n",
        sep = "")
    return(invisible(x))
  }
  for (k in seq_len(x$p)) {
    has_k <- (!is.null(x$a) && !is.null(x$a[[k]])) ||
             (!is.null(x$b) && !is.null(x$b[[k]])) ||
             (!is.null(x$W) && !is.null(x$W[[k]]))
    if (!has_k) next
    cat(sprintf("\n  coord k = %d:\n", k))
    if (!is.null(x$a) && !is.null(x$a[[k]]) && nrow(x$a[[k]]) > 0L) {
      cat("    a:\n")
      df_a <- if (level == "full") x$a[[k]] else
        x$a[[k]][, c("term", "mean"), drop = FALSE]
      print(format_coef_df(df_a, digits = digits), row.names = FALSE)
    }
    if (!is.null(x$b) && !is.null(x$b[[k]]) && nrow(x$b[[k]]) > 0L) {
      cat("    b:\n")
      df_b <- if (level == "full") x$b[[k]] else
        x$b[[k]][, c("term", "mean"), drop = FALSE]
      print(format_coef_df(df_b, digits = digits), row.names = FALSE)
    }
    if (!is.null(x$W) && !is.null(x$W[[k]]) && nrow(x$W[[k]]) > 0L) {
      cat("    W:\n")
      df_W <- if (level == "full") x$W[[k]] else
        x$W[[k]][, c("basis_idx", "x_name", "mean"), drop = FALSE]
      print(format_coef_df(df_W, digits = digits), row.names = FALSE)
    }
  }
  invisible(x)
}

#' Internal: format numeric columns of a coef data.frame for printing
#'
#' @keywords internal
#' @noRd
format_coef_df <- function(df, digits) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  num_cols <- vapply(df, is.numeric, logical(1L))
  for (col in names(df)[num_cols]) {
    if (col %in% c("g", "k", "basis_idx")) next
    df[[col]] <- formatC(df[[col]], digits = digits, format = "g",
                         flag = "-")
  }
  df
}

#' Summary method for gdpar_coef objects
#'
#' Returns a compact list of aggregated statistics: number of
#' coordinates, count of active components per type, and average of the
#' posterior means across coordinates for theta_ref.
#'
#' @param object Object of class \code{gdpar_coef}.
#' @param ... Unused.
#' @return A list with elements \code{p}, \code{n_active}
#'   (named integer with components a/b/W), \code{theta_ref_mean}
#'   (mean of theta_ref posterior means across coords),
#'   \code{summary_stats}.
#' @export
summary.gdpar_coef <- function(object, ...) {
  list(
    p              = object$p,
    n_active       = c(
      a = count_active_coords(object$a),
      b = count_active_coords(object$b),
      W = count_active_coords(object$W)
    ),
    theta_ref_mean = mean(object$theta_ref$mean),
    summary_stats  = object$summary_stats
  )
}

#' Coerce a gdpar_coef object to a long-tidy data.frame
#'
#' Flattens the hierarchical per-component, per-coordinate structure
#' into a single data.frame with columns
#' \code{(component, k, identifier, x_name, mean, q05, q50, q95)}.
#' The \code{identifier} column carries the \code{term} for \code{a}
#' and \code{b} slots, the \code{basis_idx} (formatted as a string) for
#' \code{W}, and \code{NA} for \code{theta_ref}. The \code{x_name}
#' column is \code{NA} except for \code{W} rows. Useful as input to
#' \pkg{dplyr} / \pkg{ggplot2} pipelines.
#'
#' @param x Object of class \code{gdpar_coef}.
#' @param row.names Ignored; required by the generic.
#' @param optional Ignored; required by the generic.
#' @param ... Unused.
#' @return A data.frame with one row per scalar coefficient summarized.
#' @export
as.data.frame.gdpar_coef <- function(x, row.names = NULL,
                                     optional = FALSE, ...) {
  rows <- list()
  has_group <- "g" %in% names(x$theta_ref)
  for (i in seq_len(nrow(x$theta_ref))) {
    g_val <- if (has_group) as.integer(x$theta_ref$g[i]) else NA_integer_
    rows[[length(rows) + 1L]] <- data.frame(
      component  = "theta_ref",
      g          = g_val,
      k          = x$theta_ref$k[i],
      identifier = NA_character_,
      x_name     = NA_character_,
      mean       = x$theta_ref$mean[i],
      q05        = x$theta_ref$q05[i],
      q50        = x$theta_ref$q50[i],
      q95        = x$theta_ref$q95[i],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  add_hyper <- function(slot, comp_name) {
    if (is.null(slot)) return()
    for (i in seq_len(nrow(slot))) {
      rows[[length(rows) + 1L]] <<- data.frame(
        component  = comp_name,
        g          = NA_integer_,
        k          = as.integer(slot$k[i]),
        identifier = NA_character_,
        x_name     = NA_character_,
        mean       = slot$mean[i],
        q05        = slot$q05[i],
        q50        = slot$q50[i],
        q95        = slot$q95[i],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  add_hyper(x$mu_theta_ref, "mu_theta_ref")
  add_hyper(x$sigma_theta_ref, "sigma_theta_ref")
  add_terms <- function(slot, comp_name) {
    if (is.null(slot)) return()
    for (k in seq_along(slot)) {
      df_k <- slot[[k]]
      if (is.null(df_k) || nrow(df_k) == 0L) next
      for (i in seq_len(nrow(df_k))) {
        rows[[length(rows) + 1L]] <<- data.frame(
          component  = comp_name,
          g          = NA_integer_,
          k          = as.integer(k),
          identifier = as.character(df_k$term[i]),
          x_name     = NA_character_,
          mean       = df_k$mean[i],
          q05        = df_k$q05[i],
          q50        = df_k$q50[i],
          q95        = df_k$q95[i],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    }
  }
  add_terms(x$a, "a")
  add_terms(x$b, "b")
  if (!is.null(x$W)) {
    for (k in seq_along(x$W)) {
      df_k <- x$W[[k]]
      if (is.null(df_k) || nrow(df_k) == 0L) next
      for (i in seq_len(nrow(df_k))) {
        rows[[length(rows) + 1L]] <- data.frame(
          component  = "W",
          g          = NA_integer_,
          k          = as.integer(k),
          identifier = as.character(df_k$basis_idx[i]),
          x_name     = as.character(df_k$x_name[i]),
          mean       = df_k$mean[i],
          q05        = df_k$q05[i],
          q50        = df_k$q50[i],
          q95        = df_k$q95[i],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      component  = character(0),
      g          = integer(0),
      k          = integer(0),
      identifier = character(0),
      x_name     = character(0),
      mean       = numeric(0),
      q05        = numeric(0),
      q50        = numeric(0),
      q95        = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

#' One-line formatter for gdpar_coef objects
#'
#' Useful for logs and condensed printouts.
#'
#' @param x Object of class \code{gdpar_coef}.
#' @param ... Unused.
#' @return A length-1 character vector.
#' @export
format.gdpar_coef <- function(x, ...) {
  sprintf(
    "<gdpar_coef> p=%d, components=[theta_ref, a(%d/%d), b(%d/%d), W(%d/%d)]",
    x$p,
    count_active_coords(x$a), x$p,
    count_active_coords(x$b), x$p,
    count_active_coords(x$W), x$p
  )
}

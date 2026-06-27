#' Compare an Empirical-Bayes fit against a Fully-Bayes fit
#'
#' Sub-phase 8.6.E (Charter Section 3.5, decision 2.5 Trio of vignettes).
#' Reports the operational comparison promised in v07 Section 11 between
#' a \code{gdpar_eb_fit} (Empirical Bayes via \code{\link{gdpar_eb}})
#' and a \code{gdpar_fit} (Fully Bayes via \code{\link{gdpar}}) fitted
#' on the same dataset: per-component differences in the population
#' anchor \eqn{\theta_{\text{ref}}}, the empirical total variation
#' distance between the lower-level posteriors of
#' \eqn{\xi = (a, b, W, \text{dispersion})} marginally per parameter,
#' and the operational verification of the higher-order coverage
#' discrepancy of v07 Section 6 (Proposition 7B scalar / 7B* matricial
#' / 7B* tensorial) on the nominal EB and FB credible intervals.
#'
#' The comparator is descriptive: it does not assert algorithmic
#' equivalence, nor does it test hypotheses across the EB and FB
#' inferential frames. The TV distance is computed marginally
#' parameter by parameter via histogram-based plug-in (relative bin
#' counts at a common breakpoint grid); a finite-sample correction is
#' not applied. Joint TV across the high-dimensional \eqn{\xi}
#' typically requires kernel Stein discrepancy or similar
#' density-free metrics that are out of scope of the initial 8.6.E
#' iteration; the marginal TV reported here is the operational proxy
#' recommended in v07 Section 11.1.
#'
#' @param eb_fit An object of class \code{gdpar_eb_fit} produced by
#'   \code{\link{gdpar_eb}}. Covers all four path regimes (K = 1 + p
#'   = 1; Path A K = 1 + p > 1; Path B K > 1 + p = 1; Path C K > 1 +
#'   p > 1 via the K x p tensor extension of Sub-phase 8.6.D).
#' @param fb_fit An object of class \code{gdpar_fit} produced by
#'   \code{\link{gdpar}}. Must have been fitted on the same dataset
#'   as \code{eb_fit} (same outcome, same covariates, same K / p
#'   regime). The comparator does not refit either model.
#' @param level Numeric scalar in (0, 1); credible-interval level for
#'   the coverage discrepancy reporting. Defaults to 0.95.
#' @param tv_bins Integer scalar; number of histogram bins used to
#'   approximate the marginal TV distance per parameter. Defaults to
#'   30. Larger values give a finer empirical TV but require more
#'   draws per parameter for stability.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return An object of class \code{gdpar_eb_fb_comparison} with
#'   components \code{theta_diff_table}, \code{tv_table},
#'   \code{coverage_table}, \code{level}, \code{tv_bins},
#'   \code{n_common_params}, \code{path_eb}, \code{path_fb},
#'   \code{call}, \code{warnings} (character vector recording
#'   per-helper fallback notifications: silent extraction failures of
#'   the FB \code{theta_ref} draws, missing EB or FB \eqn{\xi} draws,
#'   or zero-common-parameter TV inputs; empty in the happy path), and
#'   \code{meta}. See \code{\link{print.gdpar_eb_fb_comparison}} and
#'   \code{\link{summary.gdpar_eb_fb_comparison}} for the companion
#'   S3 methods.
#'
#' @section Path coverage:
#' The comparator handles all four EB regimes uniformly by extracting
#' the per-element anchor estimate (vector / matrix / 3D array) and
#' the corresponding lower-level posterior draws via the canonical
#' \code{posterior::as_draws_matrix} interface. For Path C
#' (K > 1 + p > 1) the \code{theta_ref_kp_hat} tensor is flattened to
#' a length-K*p vector keyed by (slot, coord) for per-element
#' comparison; the joint K x p inflation tensor is reported in the
#' coverage_table per cell via its diagonal block entries.
#'
#' @section Diagnostic value:
#' Under the standing hypotheses of v07 Section 4 (EB-MARG-ID +
#' PRIOR-FB-WEAK + HIER-COMPLEX), Theorem 7A predicts marginal TV ->
#' 0 in probability as n -> Inf (specializing to Theorems 7A* /
#' 7C* / 7C* compound multi-slot under Path A / Path B / Path C of
#' v07b Sections 4-6). The empirical TV reported here is the
#' operational diagnostic of this theoretical prediction. Persistent
#' large marginal TV across \eqn{\xi} suggests one of the discrepancy
#' conditions of Proposition 7D (multi-modality of the marginal
#' likelihood, near-singular Fisher information, informative prior on
#' theta_ref, deep hierarchy). The coverage_table operationally
#' verifies the \eqn{O(n^{-1})} under-cover claim of Proposition 7B by
#' comparing EB-nominal vs FB-nominal IC widths per anchor cell.
#'
#' @references
#' Petrone, S., Rousseau, J., and Scricciolo, C. (2014). Bayes and
#' empirical Bayes: do they merge? \emph{Biometrika}, 101(2), 285--302.
#'
#' Rousseau, J., and Szabo, B. (2017). Asymptotic behaviour of the
#' empirical Bayes posteriors associated to maximum marginal
#' likelihood estimator. \emph{Annals of Statistics}, 45(2), 833--865.
#'
#' Carlin, B. P., and Gelfand, A. E. (1990). Approaches for empirical
#' Bayes confidence intervals. \emph{JASA}, 85(409), 105--114.
#'
#' @seealso \code{\link{gdpar_eb}}, \code{\link{gdpar}},
#'   \code{vignette("v07_eb_vs_fb", package = "gdpar")},
#'   \code{vignette("v07b_eb_multivariate", package = "gdpar")},
#'   \code{vignette("vop07_eb_workflow", package = "gdpar")}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     requireNamespace("posterior", quietly = TRUE)) {
#'   set.seed(20260526L)
#'   n <- 120L
#'   df <- data.frame(x = stats::rnorm(n))
#'   df$y <- 0.5 + 0.4 * df$x + stats::rnorm(n, sd = 0.3)
#'   spec <- amm_spec(a = ~ x)
#'   fit_eb <- gdpar_eb(
#'     formula = y ~ x, family = gdpar_family("gaussian"),
#'     amm = spec, data = df,
#'     iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
#'     refresh = 0L, verbose = FALSE, seed = 1L
#'   )
#'   fit_fb <- gdpar(
#'     formula = y ~ x, family = gdpar_family("gaussian"),
#'     amm = spec, data = df,
#'     iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
#'     refresh = 0L, verbose = FALSE, seed = 1L
#'   )
#'   cmp <- gdpar_compare_eb_fb(fit_eb, fit_fb)
#'   print(cmp)
#' }
#' }
#'
#' @export
gdpar_compare_eb_fb <- function(eb_fit, fb_fit, level = 0.95,
                                 tv_bins = 30L, ...) {
  assert_inherits(eb_fit, "gdpar_eb_fit", "eb_fit")
  assert_inherits(fb_fit, "gdpar_fit", "fb_fit")
  if (!is.numeric(level) || length(level) != 1L ||
      level <= 0 || level >= 1) {
    gdpar_abort(
      "Argument 'level' must be a single numeric value in (0, 1).",
      class = "gdpar_input_error"
    )
  }
  if (!is.numeric(tv_bins) || length(tv_bins) != 1L ||
      tv_bins < 5L) {
    gdpar_abort(
      "Argument 'tv_bins' must be an integer scalar >= 5.",
      class = "gdpar_input_error"
    )
  }
  tv_bins <- as.integer(tv_bins)
  require_suggested("posterior",
                    "extract posterior draws for the EB-vs-FB comparison")
  call <- match.call()
  path_eb <- eb_fit$path %||% "eb"
  path_fb <- fb_fit$path %||% "fb"

  # D51 (Sesion 18 Etapa 2 of Sesion 8.4): accumulate fallback
  # notifications from the three helper extractors. Each helper has a
  # silent NULL fallback by design (downstream tables degrade
  # gracefully to NA), but D48 canonizes the obligation to notify the
  # user when a fallback is used and the reason it was used. The
  # orchestrator inspects the helpers' outputs and emits structured
  # gdpar_warn invocations plus persists the messages into
  # \code{$warnings} of the returned object.
  warnings_msg <- character(0L)
  emit <- function(msg, cls = "gdpar_diagnostic_warning") {
    warnings_msg <<- c(warnings_msg, msg)
    gdpar_warn(msg, class = cls)
  }

  theta_diff <- .gdpar_eb_fb_theta_diff_table(eb_fit, fb_fit, level)
  # Detect FB theta_ref-draws extraction silent failure (manifests as
  # all-NA `fb_mean` column when the FB fit either lacks draws or has a
  # template variable convention the extractor does not recognize).
  if (!is.null(theta_diff) && "fb_mean" %in% names(theta_diff) &&
      all(is.na(theta_diff$fb_mean))) {
    emit(paste(
      "Empirical-Bayes vs Fully-Bayes anchor diff table: every FB cell",
      "is NA. The extractor could not recover theta_ref draws from the",
      "FB fit (unknown template convention or empty draws). Falling",
      "back to NA-filled fb_mean / fb_se / diff columns; the EB-side",
      "estimates remain intact. See $theta_diff_table."
    ))
  }

  draws_eb <- .gdpar_eb_fb_extract_xi_draws(eb_fit$conditional_fit)
  if (is.null(draws_eb)) {
    emit(paste(
      "Empirical-Bayes xi draws unavailable from",
      "eb_fit$conditional_fit. Falling back to NULL tv_table; the",
      "anchor diff and coverage tables remain available."
    ))
  }
  draws_fb <- .gdpar_eb_fb_extract_xi_draws(fb_fit$conditional_fit %||%
                                              fb_fit$fit)
  if (is.null(draws_fb)) {
    emit(paste(
      "Fully-Bayes xi draws unavailable from fb_fit$conditional_fit /",
      "fb_fit$fit. Falling back to NULL tv_table; the anchor diff and",
      "coverage tables remain available."
    ))
  }
  tv_tbl <- .gdpar_eb_fb_tv_table(draws_eb, draws_fb, tv_bins)
  if (is.null(tv_tbl) && !is.null(draws_eb) && !is.null(draws_fb)) {
    emit(paste(
      "Marginal TV table: zero common parameter names between EB and",
      "FB xi draws. Falling back to NULL tv_table. This typically",
      "signals divergent template variable naming conventions between",
      "the two fits."
    ))
  }
  coverage_tbl <- .gdpar_eb_fb_coverage_table(eb_fit, fb_fit, level)
  if (!is.null(coverage_tbl) && "fb_width" %in% names(coverage_tbl) &&
      all(is.na(coverage_tbl$fb_width))) {
    emit(paste(
      "Coverage discrepancy table: every FB cell is NA. The extractor",
      "could not recover theta_ref draws from the FB fit for the",
      "coverage calculation. Falling back to NA-filled fb_width /",
      "width_ratio columns; the EB-side widths remain intact."
    ))
  }
  out <- list(
    theta_diff_table  = theta_diff,
    tv_table          = tv_tbl,
    coverage_table    = coverage_tbl,
    level             = level,
    tv_bins           = tv_bins,
    n_common_params   = if (is.null(tv_tbl)) 0L else nrow(tv_tbl),
    path_eb           = path_eb,
    path_fb           = path_fb,
    family_eb         = eb_fit$family$name %||% NA_character_,
    family_fb         = fb_fit$family$name %||% NA_character_,
    call              = call,
    warnings          = warnings_msg,
    meta = list(
      mode = "compare_eb_fb",
      note = paste(
        "Marginal TV distance per parameter via plug-in histogram",
        "(tv_bins =", tv_bins, "). Coverage discrepancy reports",
        "EB-nominal vs FB-nominal IC widths per anchor cell at",
        "level =", level, "."
      )
    )
  )
  class(out) <- c("gdpar_eb_fb_comparison", "list")
  out
}

#' Internal: build the theta_diff_table comparing EB and FB anchors
#'
#' Returns a data frame with one row per anchor cell. Path C produces
#' rows keyed by (group, slot, coord); Path B by (group, slot); Path A
#' by (group, coord); K = 1 + p = 1 by group only. The columns
#' \code{eb_estimate}, \code{eb_se}, \code{fb_mean}, \code{fb_se},
#' \code{diff} (\code{eb_estimate - fb_mean}), and \code{diff_rel}
#' (\code{diff / fb_se}) are filled for every row.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_fb_theta_diff_table <- function(eb_fit, fb_fit, level) {
  is_path_C <- identical(eb_fit$path, "eb_KxP")
  fb_draws <- tryCatch(
    .gdpar_eb_fb_extract_theta_ref_draws_fb(fb_fit),
    error = function(e) NULL
  )
  if (is_path_C) {
    J <- dim(eb_fit$theta_ref_kp_hat)[1L]
    K <- eb_fit$K
    p <- eb_fit$p
    rows <- vector("list", J * K * p)
    idx <- 1L
    for (g in seq_len(J)) {
      for (k in seq_len(K)) {
        for (c in seq_len(p)) {
          eb_est <- eb_fit$theta_ref_kp_hat[g, k, c]
          eb_se <- eb_fit$theta_ref_kp_se[g, k, c]
          fb_mean <- NA_real_
          fb_se <- NA_real_
          if (!is.null(fb_draws$kp)) {
            cell <- fb_draws$kp[[g]][[k]][, c]
            if (!is.null(cell) && length(cell) > 0L) {
              fb_mean <- mean(cell)
              fb_se <- stats::sd(cell)
            }
          }
          rows[[idx]] <- data.frame(
            group       = g,
            slot        = eb_fit$slot_names[k],
            coord       = c,
            eb_estimate = eb_est,
            eb_se       = eb_se,
            fb_mean     = fb_mean,
            fb_se       = fb_se,
            diff        = eb_est - fb_mean,
            diff_rel    = if (isTRUE(is.finite(fb_se) && fb_se > 0))
                            (eb_est - fb_mean) / fb_se else NA_real_,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
    return(do.call(rbind, rows))
  }
  est_eb <- as.numeric(eb_fit$theta_ref_hat)
  se_eb <- as.numeric(eb_fit$theta_ref_se)
  n_cells <- length(est_eb)
  fb_means <- rep(NA_real_, n_cells)
  fb_ses <- rep(NA_real_, n_cells)
  if (!is.null(fb_draws$flat)) {
    n_min <- min(n_cells, length(fb_draws$flat$means))
    fb_means[seq_len(n_min)] <- fb_draws$flat$means[seq_len(n_min)]
    fb_ses[seq_len(n_min)] <- fb_draws$flat$ses[seq_len(n_min)]
  }
  data.frame(
    cell        = seq_len(n_cells),
    eb_estimate = est_eb,
    eb_se       = se_eb,
    fb_mean     = fb_means,
    fb_se       = fb_ses,
    diff        = est_eb - fb_means,
    diff_rel    = ifelse(is.finite(fb_ses) & fb_ses > 0,
                          (est_eb - fb_means) / fb_ses, NA_real_),
    stringsAsFactors = FALSE
  )
}

#' Internal: extract theta_ref draws from an FB fit (path-aware)
#'
#' For K = 1 + p = 1 / Path A / Path B the FB fit stores draws under
#' \code{theta_ref[...]} or \code{theta_ref_k[...]} variable names.
#' For Path C the FB fit theoretically uses \code{theta_ref_kp[...]}
#' (the K x p FB template is itself a follow-on debt of 8.4 unification
#' per Charter and [[project_gdpar_deuda_8_4_unificacion_stan]]); this
#' helper attempts the canonical extraction and returns NULL if the
#' draws are not present in the recognized convention.
#'
#' @return A list with components \code{flat} (named list with
#'   \code{means} / \code{ses} vectors keyed by cell index) for non-Path-C
#'   regimes; or \code{kp} (nested list of matrices [J][K][, p]) for
#'   Path C. NULL when the extraction fails.
#' @keywords internal
#' @noRd
.gdpar_eb_fb_extract_theta_ref_draws_fb <- function(fb_fit) {
  draws_obj <- fb_fit$conditional_fit %||% fb_fit$fit
  if (is.null(draws_obj)) return(NULL)
  draws <- tryCatch(draws_obj$draws(), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  vars <- tryCatch(dimnames(draws)$variable %||% character(0L),
                    error = function(e) character(0L))
  if (length(vars) == 0L) return(NULL)
  # Path B convention
  vars_b <- grep("^theta_ref_k\\[", vars, value = TRUE)
  # Path C convention
  vars_c <- grep("^theta_ref_kp\\[", vars, value = TRUE)
  # K=1+p=1 / Path A convention
  vars_default <- grep("^theta_ref(\\[|$)", vars, value = TRUE)
  if (length(vars_c) > 0L) {
    # Path C extraction (kp), not consumed by the diff table when
    # eb_fit is K=1+p=1; only the kp branch above touches this.
    mat <- posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = vars_c)
    )
    return(list(kp = .gdpar_eb_fb_unpack_kp(mat, vars_c)))
  }
  if (length(vars_b) > 0L) {
    mat <- posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = vars_b)
    )
    means <- colMeans(mat)
    ses <- apply(mat, 2L, stats::sd)
    return(list(flat = list(means = unname(means),
                              ses = unname(ses))))
  }
  if (length(vars_default) > 0L) {
    mat <- posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = vars_default)
    )
    means <- colMeans(mat)
    ses <- apply(mat, 2L, stats::sd)
    return(list(flat = list(means = unname(means),
                              ses = unname(ses))))
  }
  NULL
}

#' Internal: unpack Path C theta_ref_kp[g,k,c] draws into nested list
#'
#' @keywords internal
#' @noRd
.gdpar_eb_fb_unpack_kp <- function(mat, vars_c) {
  parsed <- regmatches(vars_c, regexec("\\[(\\d+),(\\d+),(\\d+)\\]",
                                         vars_c))
  triples <- do.call(rbind, lapply(parsed, function(x) {
    if (length(x) < 4L) c(NA, NA, NA) else as.integer(x[2:4])
  }))
  J <- max(triples[, 1L], na.rm = TRUE)
  K <- max(triples[, 2L], na.rm = TRUE)
  p <- max(triples[, 3L], na.rm = TRUE)
  kp <- vector("list", J)
  for (g in seq_len(J)) {
    kp[[g]] <- vector("list", K)
    for (k in seq_len(K)) {
      kp[[g]][[k]] <- matrix(NA_real_, nrow = nrow(mat), ncol = p)
      for (c in seq_len(p)) {
        nm <- sprintf("theta_ref_kp[%d,%d,%d]", g, k, c)
        if (nm %in% colnames(mat)) {
          kp[[g]][[k]][, c] <- mat[, nm]
        }
      }
    }
  }
  kp
}

#' Internal: extract the xi posterior draws (a_coef / c_b / W / dispersion)
#'
#' Used by both the EB conditional fit and the FB fit. Filters out
#' \code{lp__}, \code{theta_ref*}, raw / packed / per-coord helper
#' variables, and the generated quantities (\code{eta*},
#' \code{log_lik*}, \code{y_pred*}, \code{theta_i*}). Returns a
#' \code{posterior::draws_matrix} subset on the retained variables.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_fb_extract_xi_draws <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  draws <- tryCatch(fit_obj$draws(), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  vars <- tryCatch(dimnames(draws)$variable %||% character(0L),
                    error = function(e) character(0L))
  if (length(vars) == 0L) return(NULL)
  keep <- vars[!grepl(
    paste(
      "^lp__$|^theta_ref|^mu_theta_ref|^sigma_theta_ref",
      "^eta|^eta_kp|^log_lik|^y_pred|^theta_i",
      "^a_raw|^c_b_raw|^c_b_kp_raw|^W_raw",
      sep = "|"
    ),
    vars
  )]
  if (length(keep) == 0L) return(NULL)
  posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = keep)
  )
}

#' Internal: marginal TV distance per parameter via histogram plug-in
#'
#' For each common parameter between draws_eb and draws_fb, computes
#' an empirical TV distance from the bin counts of the two marginals
#' over a shared support grid: \code{TV = 0.5 * sum(|p_eb - p_fb|)}
#' over the bins. Returns a data frame with columns
#' \code{parameter}, \code{tv}, \code{n_eb}, \code{n_fb},
#' \code{mean_eb}, \code{mean_fb}.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_fb_tv_table <- function(draws_eb, draws_fb, tv_bins) {
  if (is.null(draws_eb) || is.null(draws_fb)) return(NULL)
  common <- intersect(colnames(draws_eb), colnames(draws_fb))
  if (length(common) == 0L) return(NULL)
  rows <- vector("list", length(common))
  for (i in seq_along(common)) {
    nm <- common[i]
    x_eb <- as.numeric(draws_eb[, nm])
    x_fb <- as.numeric(draws_fb[, nm])
    rng <- range(c(x_eb, x_fb), na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) {
      tv <- NA_real_
    } else {
      breaks <- seq(rng[1L], rng[2L], length.out = tv_bins + 1L)
      h_eb <- graphics::hist(x_eb, breaks = breaks, plot = FALSE)$counts
      h_fb <- graphics::hist(x_fb, breaks = breaks, plot = FALSE)$counts
      p_eb <- h_eb / sum(h_eb)
      p_fb <- h_fb / sum(h_fb)
      tv <- 0.5 * sum(abs(p_eb - p_fb))
    }
    rows[[i]] <- data.frame(
      parameter = nm,
      tv        = tv,
      n_eb      = length(x_eb),
      n_fb      = length(x_fb),
      mean_eb   = mean(x_eb),
      mean_fb   = mean(x_fb),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Internal: build the coverage discrepancy table over theta_ref cells
#'
#' Reports per anchor cell: EB credible-interval width (with inflation
#' applied when \code{eb_fit$correction_applied}), FB credible-interval
#' width derived from the FB posterior draws, and the ratio
#' \code{width_eb / width_fb} that operationally diagnoses the
#' \eqn{O(n^{-1})} under-cover prediction of v07 Section 6 Proposition
#' 7B (Path A / Path B / Path C extensions of v07b Section 5).
#'
#' @keywords internal
#' @noRd
.gdpar_eb_fb_coverage_table <- function(eb_fit, fb_fit, level) {
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  is_path_C <- identical(eb_fit$path, "eb_KxP")
  fb_draws <- tryCatch(
    .gdpar_eb_fb_extract_theta_ref_draws_fb(fb_fit),
    error = function(e) NULL
  )
  if (is_path_C) {
    J <- dim(eb_fit$theta_ref_kp_hat)[1L]
    K <- eb_fit$K
    p <- eb_fit$p
    rows <- vector("list", J * K * p)
    idx <- 1L
    correction_on <- isTRUE(eb_fit$correction_applied)
    tensor <- eb_fit$correction_tensor_constant
    for (g in seq_len(J)) {
      for (k in seq_len(K)) {
        for (c in seq_len(p)) {
          eb_se <- eb_fit$theta_ref_kp_se[g, k, c]
          inflate_kc <- if (correction_on && !is.null(tensor) &&
                            all(is.finite(tensor[k, c, c]))) {
                          sqrt(1 + tensor[k, c, c] / max(1, J))
                        } else {
                          1
                        }
          eb_width <- 2 * z * eb_se * inflate_kc
          fb_width <- NA_real_
          if (!is.null(fb_draws$kp)) {
            cell <- fb_draws$kp[[g]][[k]][, c]
            if (!is.null(cell) && length(cell) > 0L) {
              fb_width <- 2 * z * stats::sd(cell)
            }
          }
          rows[[idx]] <- data.frame(
            group        = g,
            slot         = eb_fit$slot_names[k],
            coord        = c,
            eb_width     = eb_width,
            fb_width     = fb_width,
            width_ratio  = if (isTRUE(is.finite(fb_width) &&
                                       fb_width > 0))
                              eb_width / fb_width else NA_real_,
            inflation    = inflate_kc,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
    return(do.call(rbind, rows))
  }
  est_eb <- as.numeric(eb_fit$theta_ref_hat)
  se_eb <- as.numeric(eb_fit$theta_ref_se)
  inflate <- if (isTRUE(eb_fit$correction_applied))
               sqrt(1 + (eb_fit$eb_correction_constant %||% 0) /
                            max(1, length(est_eb)))
             else 1
  eb_widths <- 2 * z * se_eb * inflate
  fb_widths <- rep(NA_real_, length(est_eb))
  if (!is.null(fb_draws$flat)) {
    n_min <- min(length(est_eb), length(fb_draws$flat$ses))
    fb_widths[seq_len(n_min)] <- 2 * z * fb_draws$flat$ses[seq_len(n_min)]
  }
  data.frame(
    cell        = seq_along(est_eb),
    eb_width    = eb_widths,
    fb_width    = fb_widths,
    width_ratio = ifelse(is.finite(fb_widths) & fb_widths > 0,
                          eb_widths / fb_widths, NA_real_),
    inflation   = inflate,
    stringsAsFactors = FALSE
  )
}

#' Compare a fit snapshot against a golden reference (experimental)
#'
#' Four-layer comparator that locks the posterior summaries and sampler
#' diagnostics of a fitted \code{gdpar} model against a persisted
#' reference snapshot. Designed for regression testing of MCMC outputs:
#' any future deviation that exceeds the principled tolerance is
#' flagged as a failure together with the layer, item, expected and
#' observed values, and the diagnostic delta. Layer A (structural)
#' compares class signatures, slot shapes and column names; layer B
#' (discrete) enforces bit-exact agreement on integer sampler
#' diagnostics; layer C (continuous) uses Monte Carlo standard error
#' from the golden as the principled tolerance, \code{|obs - exp| <=
#' k_sigma * MC_SE_golden}; layer D (sanity) checks absolute floors
#' that must hold regardless of the golden (R-hat, ESS, divergent
#' percentage, E-BFMI). All four layers are evaluated; failures are
#' aggregated rather than short-circuited so that the caller obtains a
#' comprehensive report.
#'
#' This function is the comparator counterpart of
#' \code{\link{gdpar_snapshot_fit}}. The typical workflow is documented
#' in \code{vignette("vop03_regression_testing", package = "gdpar")}.
#'
#' @param observed A list with the snapshot schema produced by
#'   \code{\link{gdpar_snapshot_fit}}: fields \code{structural},
#'   \code{discrete}, \code{continuous}, \code{sanity}, and
#'   \code{parametrization_resolved}.
#' @param golden A list of identical schema, typically loaded via
#'   \code{readRDS()} from an archived reference snapshot.
#' @param k_sigma Numeric scalar; tolerance multiplier for layer C
#'   (continuous). Default 3. Lower values increase sensitivity to
#'   posterior mean drift; higher values relax it.
#' @param sanity_floor Optional named list overriding the default
#'   absolute thresholds for layer D. Defaults to
#'   \code{list(rhat_max = 1.05, ess_bulk_min = 100, ess_tail_min = 100,
#'   divergent_pct = 0.01, ebfmi_min = 0.3)}.
#'
#' @return A list with components \code{passed} (logical scalar;
#'   \code{TRUE} iff every layer is clean), \code{failures} (data.frame
#'   with one row per failure; columns \code{layer}, \code{item},
#'   \code{expected}, \code{observed}, \code{delta}, \code{threshold},
#'   \code{severity}), and \code{by_layer} (named integer vector of
#'   failure counts per layer).
#'
#' @section Status:
#' This function is flagged \strong{experimental}. The snapshot schema
#' is versioned at \code{schema_version = 1L}; future Blocks 6-9 of the
#' development roadmap may add fields and bump the schema. The
#' tolerance contract (\code{k_sigma}, sanity floors) is also subject
#' to refinement until the first stable release.
#'
#' @examples
#' \donttest{
#' make_snapshot <- function(mean_val) list(
#'   structural = NULL,
#'   discrete   = NULL,
#'   continuous = list(theta_ref = list(
#'     "theta_ref[1]" = list(mean = mean_val, sd = 0.01,
#'                            ess_bulk = 1000, ess_tail = 1000,
#'                            rhat = 1.001, mc_se = 0.001)
#'   )),
#'   sanity     = list(rhat_max = 1.0, ess_bulk_min = 1000,
#'                      ess_tail_min = 1000, divergent_pct = 0,
#'                      ebfmi_min = 1.0)
#' )
#' golden   <- make_snapshot(0.5)
#' observed <- make_snapshot(0.5005)
#' cmp <- gdpar_golden_compare(observed, golden, k_sigma = 3)
#' cmp$passed
#' }
#'
#' @seealso \code{\link{gdpar_snapshot_fit}}, \code{\link{gdpar}}.
#' @keywords experimental
#' @export
gdpar_golden_compare <- function(observed, golden, k_sigma = 3,
                                  sanity_floor = NULL) {
  if (!is.list(observed) || !is.list(golden)) {
    gdpar_abort(
      "gdpar_golden_compare(): 'observed' and 'golden' must both be lists.",
      class = "gdpar_internal_error"
    )
  }
  failures <- list()
  add <- function(layer, item, expected, observed_v, delta = NA_real_,
                  threshold = NA_real_, severity = "fail") {
    failures[[length(failures) + 1L]] <<- data.frame(
      layer     = layer,
      item      = item,
      expected  = format_compact(expected),
      observed  = format_compact(observed_v),
      delta     = as.numeric(delta),
      threshold = as.numeric(threshold),
      severity  = severity,
      stringsAsFactors = FALSE
    )
  }
  golden_compare_structural(observed$structural, golden$structural, add)
  golden_compare_discrete(observed$discrete, golden$discrete, add)
  golden_compare_continuous(observed$continuous, golden$continuous,
                             k_sigma, add)
  defaults <- list(rhat_max = 1.05, ess_bulk_min = 100,
                    ess_tail_min = 100, divergent_pct = 0.01,
                    ebfmi_min = 0.3)
  if (is.null(sanity_floor)) sanity_floor <- defaults
  else sanity_floor <- utils::modifyList(defaults, sanity_floor)
  golden_compare_sanity(observed$sanity, sanity_floor, add)

  if (length(failures) == 0L) {
    failures_df <- data.frame(
      layer = character(0), item = character(0),
      expected = character(0), observed = character(0),
      delta = numeric(0), threshold = numeric(0),
      severity = character(0), stringsAsFactors = FALSE
    )
  } else {
    failures_df <- do.call(rbind, failures)
  }
  by_layer <- table(failures_df$layer)
  list(
    passed = nrow(failures_df) == 0L,
    failures = failures_df,
    by_layer = as.integer(by_layer)[
      match(c("structural", "discrete", "continuous", "sanity"),
            names(by_layer))
    ] |> stats::setNames(c("structural", "discrete", "continuous",
                            "sanity"))
  )
}

#' Internal: format any value to compact one-line string for the report
#'
#' @keywords internal
#' @noRd
format_compact <- function(x) {
  if (is.null(x)) return("<NULL>")
  if (length(x) == 0L) return("<empty>")
  if (is.list(x)) {
    return(paste0("<list:", length(x), ">"))
  }
  if (is.character(x)) {
    return(paste(x, collapse = ","))
  }
  if (is.logical(x)) {
    return(paste(as.character(x), collapse = ","))
  }
  if (is.numeric(x)) {
    if (length(x) <= 4L) {
      return(paste(trimws(formatC(x, digits = 6, format = "g")),
                    collapse = ","))
    } else {
      return(sprintf("<num[%d]>", length(x)))
    }
  }
  paste(as.character(x), collapse = ",")
}

#' Build a reproducibility snapshot of a fitted gdpar model (experimental)
#'
#' Extract the four-layer snapshot consumed by
#' \code{\link{gdpar_golden_compare}}: posterior summaries with Monte
#' Carlo standard error per variable, integer sampler diagnostics,
#' aggregated sanity floors, and the structural class signature of the
#' fit. The snapshot is the canonical reference object for regression
#' testing of MCMC outputs across package upgrades, cmdstan upgrades,
#' or refactors that touch the sampling side of the package.
#'
#' Persist the returned list via \code{saveRDS()} to lock the fit, and
#' compare against future fits via \code{\link{gdpar_golden_compare}}.
#'
#' @param fit A \code{gdpar_fit} object returned by \code{\link{gdpar}}.
#'
#' @return A list with fields \code{structural} (class signatures, slot
#'   shapes, column names), \code{discrete} (integer sampler
#'   diagnostics: \code{n_divergent}, \code{treedepth_max_n},
#'   \code{treedepth_max_value},
#'   \code{n_leapfrog_total_per_chain}, \code{ebfmi_min}),
#'   \code{continuous} (per-variable posterior mean, sd, ESS bulk /
#'   tail, R-hat, Monte Carlo standard error), \code{sanity}
#'   (aggregated convergence floors), and
#'   \code{parametrization_resolved} (resolved CP/NCP flags and the
#'   aggregation method used by the pre-flight).
#'
#' @section Status:
#' This function is flagged \strong{experimental}. The schema is at
#' \code{schema_version = 1L}; future Blocks 6-9 of the development
#' roadmap may add fields and bump the schema with documented
#' migration.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   set.seed(1L)
#'   df <- data.frame(x1 = rnorm(50), y = rnorm(50))
#'   fit <- gdpar(
#'     y ~ x1, amm = amm_spec(a = ~ x1), data = df,
#'     chains = 1L, iter_warmup = 100L, iter_sampling = 100L,
#'     refresh = 0L, verbose = FALSE, seed = 1L
#'   )
#'   snap <- gdpar_snapshot_fit(fit)
#'   names(snap)
#' }
#' }
#'
#' @seealso \code{\link{gdpar_golden_compare}}, \code{\link{gdpar}}.
#' @keywords experimental
#' @export
gdpar_snapshot_fit <- function(fit) {
  assert_inherits(fit, "gdpar_fit", "fit")
  draws_arr <- fit$fit$draws(format = "draws_array")
  continuous <- list()
  for (vn in c("theta_ref", "a_coef", "sigma_a", "sigma_y", "lp__")) {
    s <- summarize_draws_var_golden(draws_arr, vn)
    if (!is.null(s)) continuous[[vn]] <- s
  }
  discrete <- extract_discrete_diag_golden(fit)
  sanity <- extract_sanity_golden(draws_arr, discrete)
  structural <- build_structural_snapshot_golden(fit)
  list(
    structural               = structural,
    discrete                 = discrete,
    continuous               = continuous,
    sanity                   = sanity,
    parametrization_resolved = list(
      cp_a        = fit$parametrization$cp_a,
      cp_W        = fit$parametrization$cp_W,
      cp_a_per_k  = fit$parametrization$cp_a_per_k,
      cp_W_per_k  = fit$parametrization$cp_W_per_k,
      aggregation = fit$parametrization$meta$aggregation
    )
  )
}

#' Internal: per-variable posterior summary for the golden snapshot
#'
#' Computes per-element mean, sd, ess_bulk, ess_tail, rhat, and the
#' Monte Carlo standard error (\code{sd / sqrt(ess_bulk)}) for every
#' element of a variable matching the prefix \code{var_name}. Returns
#' \code{NULL} if the variable is absent.
#'
#' @keywords internal
#' @noRd
summarize_draws_var_golden <- function(draws_array, var_name) {
  vars <- posterior::variables(draws_array)
  idx <- grep(sprintf("^%s($|\\[)", var_name), vars, value = FALSE)
  if (length(idx) == 0L) return(NULL)
  sub_names <- vars[idx]
  out <- list()
  for (vn in sub_names) {
    s <- posterior::summarise_draws(
      posterior::subset_draws(draws_array, variable = vn),
      mean = ~ mean(.),
      sd = ~ stats::sd(.),
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail,
      rhat = posterior::rhat
    )
    out[[vn]] <- list(
      mean     = as.numeric(s$mean),
      sd       = as.numeric(s$sd),
      ess_bulk = as.numeric(s$ess_bulk),
      ess_tail = as.numeric(s$ess_tail),
      rhat     = as.numeric(s$rhat),
      mc_se    = as.numeric(s$sd) / sqrt(as.numeric(s$ess_bulk))
    )
  }
  out
}

#' Internal: discrete sampler diagnostics for the golden snapshot
#'
#' Extracts integer-valued sampler metrics from a \code{gdpar_fit}:
#' total divergences, count of transitions that hit
#' \code{max_treedepth}, maximum treedepth across draws, total leapfrog
#' steps per chain, and the minimum E-BFMI across chains.
#'
#' @keywords internal
#' @noRd
extract_discrete_diag_golden <- function(fit) {
  cs_fit <- fit$fit
  sampler_diag <- cs_fit$sampler_diagnostics(format = "draws_df")
  ddiag <- cs_fit$diagnostic_summary(quiet = TRUE)
  list(
    n_divergent     = as.integer(sum(ddiag$num_divergent)),
    treedepth_max_n = as.integer(sum(ddiag$num_max_treedepth)),
    n_leapfrog_total_per_chain = vapply(
      seq_along(ddiag$num_divergent),
      function(ch) {
        rows <- sampler_diag$.chain == ch
        as.integer(sum(sampler_diag$n_leapfrog__[rows]))
      },
      integer(1L)
    ),
    treedepth_max_value = as.integer(max(sampler_diag$treedepth__)),
    ebfmi_min = as.numeric(min(ddiag$ebfmi))
  )
}

#' Internal: sanity floor snapshot
#'
#' Aggregates posterior-side convergence floors (max rhat, min
#' ess_bulk, min ess_tail) and sampler-side relative floors
#' (divergent_pct, ebfmi_min). Used both for persistence and for
#' layer D of \code{gdpar_golden_compare()}.
#'
#' @keywords internal
#' @noRd
extract_sanity_golden <- function(draws_array, discrete) {
  s <- posterior::summarise_draws(
    draws_array,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  s <- s[!is.na(s$rhat), ]
  n_chains <- length(discrete$n_leapfrog_total_per_chain)
  n_sampling <- posterior::ndraws(draws_array) / n_chains
  total_transitions <- n_chains * as.integer(n_sampling)
  list(
    rhat_max      = as.numeric(max(s$rhat, na.rm = TRUE)),
    ess_bulk_min  = as.numeric(min(s$ess_bulk, na.rm = TRUE)),
    ess_tail_min  = as.numeric(min(s$ess_tail, na.rm = TRUE)),
    divergent_pct = discrete$n_divergent / total_transitions,
    ebfmi_min     = discrete$ebfmi_min
  )
}

#' Internal: structural snapshot
#'
#' Records class signatures, slot shapes, and column names so layer
#' A of \code{gdpar_golden_compare()} can detect silent API changes
#' before any numerical comparison runs.
#'
#' @keywords internal
#' @noRd
build_structural_snapshot_golden <- function(fit) {
  cf <- coef(fit)
  list(
    fit_class       = class(fit),
    coef_class      = class(cf),
    p               = as.integer(cf$p),
    summary_stats   = cf$summary_stats,
    components      = names(cf)[!names(cf) %in%
                                   c("p", "summary_stats")],
    theta_ref_cols  = colnames(cf$theta_ref),
    theta_ref_nrow  = nrow(cf$theta_ref),
    a_class         = if (is.null(cf$a)) "NULL" else class(cf$a),
    a_length        = if (is.null(cf$a)) NA_integer_
                       else length(cf$a),
    a_per_k_cols    = if (is.null(cf$a)) NULL
                       else lapply(cf$a, function(z)
                                     if (is.null(z)) NULL
                                     else colnames(z)),
    b_class         = if (is.null(cf$b)) "NULL" else class(cf$b),
    b_length        = if (is.null(cf$b)) NA_integer_
                       else length(cf$b),
    W_class         = if (is.null(cf$W)) "NULL" else class(cf$W),
    W_length        = if (is.null(cf$W)) NA_integer_
                       else length(cf$W),
    parametrization_keys = sort(names(fit$parametrization),
                                  method = "radix")
  )
}

#' Internal: layer A (structural)
#'
#' @keywords internal
#' @noRd
golden_compare_structural <- function(obs, exp, add) {
  if (is.null(obs) && is.null(exp)) return(invisible(NULL))
  if (is.null(obs) || is.null(exp)) {
    add("structural", "snapshot_present",
        !is.null(exp), !is.null(obs))
    return(invisible(NULL))
  }
  scalar_keys <- c("p", "a_length", "b_length", "W_length",
                    "theta_ref_nrow", "a_class", "b_class",
                    "W_class")
  for (key in scalar_keys) {
    if (!identical(obs[[key]], exp[[key]])) {
      add("structural", key, exp[[key]], obs[[key]])
    }
  }
  ordered_vec_keys <- c("fit_class", "coef_class", "summary_stats",
                         "theta_ref_cols")
  for (key in ordered_vec_keys) {
    if (!identical(obs[[key]], exp[[key]])) {
      add("structural", key, exp[[key]], obs[[key]])
    }
  }
  set_keys <- c("components", "parametrization_keys")
  for (key in set_keys) {
    if (!setequal(obs[[key]], exp[[key]])) {
      add("structural", key, exp[[key]], obs[[key]])
    }
  }
  if (!identical(obs$a_per_k_cols, exp$a_per_k_cols)) {
    add("structural", "a_per_k_cols",
        format_compact(exp$a_per_k_cols),
        format_compact(obs$a_per_k_cols))
  }
}

#' Internal: layer B (discrete, bit-exact integers)
#'
#' @keywords internal
#' @noRd
golden_compare_discrete <- function(obs, exp, add) {
  if (is.null(obs) && is.null(exp)) return(invisible(NULL))
  if (is.null(obs) || is.null(exp)) {
    add("discrete", "snapshot_present",
        !is.null(exp), !is.null(obs))
    return(invisible(NULL))
  }
  int_keys <- c("n_divergent", "treedepth_max_n",
                 "treedepth_max_value")
  for (key in int_keys) {
    if (!identical(as.integer(obs[[key]]), as.integer(exp[[key]]))) {
      add("discrete", key, exp[[key]], obs[[key]],
          delta = obs[[key]] - exp[[key]])
    }
  }
  if (!is.null(exp$n_leapfrog_total_per_chain) &&
        !is.null(obs$n_leapfrog_total_per_chain)) {
    if (!identical(as.integer(obs$n_leapfrog_total_per_chain),
                    as.integer(exp$n_leapfrog_total_per_chain))) {
      add("discrete", "n_leapfrog_total_per_chain",
          exp$n_leapfrog_total_per_chain,
          obs$n_leapfrog_total_per_chain)
    }
  }
}

#' Internal: layer C (continuous, k * MC_SE tolerance)
#'
#' For every parameter present in both \code{exp} and \code{obs}
#' continuous summaries, computes
#' \code{|obs$mean - exp$mean|} and flags it as a failure if it exceeds
#' \code{k_sigma * exp$mc_se}. The MC standard error is the principled
#' bound for "indistinguishable under sampling noise of the golden."
#'
#' @keywords internal
#' @noRd
golden_compare_continuous <- function(obs, exp, k_sigma, add) {
  if (is.null(obs) && is.null(exp)) return(invisible(NULL))
  if (is.null(obs) || is.null(exp)) {
    add("continuous", "snapshot_present",
        !is.null(exp), !is.null(obs))
    return(invisible(NULL))
  }
  common_groups <- intersect(names(obs), names(exp))
  for (grp in common_groups) {
    obs_g <- obs[[grp]]
    exp_g <- exp[[grp]]
    common_vars <- intersect(names(obs_g), names(exp_g))
    for (vn in common_vars) {
      o <- obs_g[[vn]]
      e <- exp_g[[vn]]
      if (is.null(o$mean) || is.null(e$mean) || is.null(e$mc_se)) {
        next
      }
      delta <- abs(o$mean - e$mean)
      thr <- k_sigma * pmax(e$mc_se, .Machine$double.eps)
      if (any(delta > thr, na.rm = TRUE)) {
        add("continuous", paste0(grp, ":", vn),
            e$mean, o$mean,
            delta = max(delta), threshold = max(thr),
            severity = "fail")
      }
    }
  }
}

#' Internal: layer D (sanity, absolute floors)
#'
#' @keywords internal
#' @noRd
golden_compare_sanity <- function(obs, floor, add) {
  if (is.null(obs)) {
    add("sanity", "snapshot_present", "non-null", "NULL")
    return(invisible(NULL))
  }
  if (!is.null(obs$rhat_max) && obs$rhat_max > floor$rhat_max) {
    add("sanity", "rhat_max", floor$rhat_max, obs$rhat_max,
        delta = obs$rhat_max - floor$rhat_max,
        threshold = floor$rhat_max)
  }
  if (!is.null(obs$ess_bulk_min) &&
        obs$ess_bulk_min < floor$ess_bulk_min) {
    add("sanity", "ess_bulk_min", floor$ess_bulk_min,
        obs$ess_bulk_min,
        delta = obs$ess_bulk_min - floor$ess_bulk_min,
        threshold = floor$ess_bulk_min)
  }
  if (!is.null(obs$ess_tail_min) &&
        obs$ess_tail_min < floor$ess_tail_min) {
    add("sanity", "ess_tail_min", floor$ess_tail_min,
        obs$ess_tail_min,
        delta = obs$ess_tail_min - floor$ess_tail_min,
        threshold = floor$ess_tail_min)
  }
  if (!is.null(obs$divergent_pct) &&
        obs$divergent_pct > floor$divergent_pct) {
    add("sanity", "divergent_pct", floor$divergent_pct,
        obs$divergent_pct,
        delta = obs$divergent_pct - floor$divergent_pct,
        threshold = floor$divergent_pct)
  }
  if (!is.null(obs$ebfmi_min) && obs$ebfmi_min < floor$ebfmi_min) {
    add("sanity", "ebfmi_min", floor$ebfmi_min, obs$ebfmi_min,
        delta = obs$ebfmi_min - floor$ebfmi_min,
        threshold = floor$ebfmi_min)
  }
}

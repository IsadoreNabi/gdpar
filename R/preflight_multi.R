#' Internal: multivariate (p > 1) pre-flight diagnostic for CP/NCP
#'
#' Counterpart of \code{preflight_parametrization()} for the
#' coord-wise factorized multivariate path. Compiles and samples a
#' short NCP version of the multivariate Stan template and applies
#' Path B' (info-ratio z-test with chain-aware block bootstrap) to
#' each coordinate of each effective component independently.
#'
#' Filter ordering mirrors the univariate path:
#' \enumerate{
#'   \item Divergence attribution per coordinate of each component.
#'         For component \code{a}, the per-\eqn{k} score uses
#'         \code{log(sigma_a[k])} and \code{||a_raw[chunk_k]||}; for
#'         component \code{W} it uses \code{log(sigma_W[1])} and the
#'         norm of the \code{k}-th block of \code{W_raw}.
#'   \item E-BFMI: global, not attributable per-component; if
#'         \code{ebfmi_min < 0.3}, every still-undecided coordinate
#'         of every effective component flips to CP.
#'   \item Info ratio (Path B'). For component \code{a},
#'         per-\eqn{k} the effective coefficient matrix is
#'         \code{a_coef[, k, 1:J_a_free[k]]} and the reference
#'         scale is \code{sigma_a[, k]}. For component \code{W},
#'         per-\eqn{k} the effective coefficient at coordinate
#'         \eqn{jj} per draw \eqn{t} is \deqn{ \mathrm{eff}[t, jj]
#'         = \sum_{m=1}^{W_{per_k_dim}} (\theta_{ref}[t, k]^m -
#'         \theta_{anchor}[k]^m) W_{raw}[t, m_k, jj] \sigma_W[t] }
#'         where \eqn{m_k = (k - 1) W_{per_k_dim} + m}; the
#'         reference scale per draw is \eqn{\sigma_W[t]
#'         \sqrt{\sum_m (\theta_{ref}[t, k]^m -
#'         \theta_{anchor}[k]^m)^2}}.
#' }
#'
#' The per-coordinate decisions are then aggregated to a
#' per-component decision via \code{aggregation}: \code{"any_ncp"}
#' (default, conservative), \code{"majority"}, or \code{"per_k"}.
#' The full per-coord, per-component report is returned as a
#' \code{gdpar_preflight_report} S3 object.
#'
#' Phase H.1 contract: the aggregation strategies \code{"any_ncp"}
#' and \code{"majority"} are wired through to the existing
#' multivariate Stan template, which honors a single (CP, NCP) pair
#' per component. The \code{"per_k"} strategy requires the per-k
#' Stan template wiring delivered in Phase H.2 and is rejected at
#' the resolver level until then.
#'
#' @param prior A \code{gdpar_prior} object.
#' @param stan_data List assembled by \code{.assemble_stan_data_multi}.
#' @param amm An \code{amm_spec} object with \code{amm$p > 1L}.
#' @param preflight_seed Integer or \code{NULL}.
#' @param verbose Logical. If \code{FALSE}, pre-flight stderr is sunk
#'   to a tempfile and discarded.
#' @param tau_cp Numeric scalar. Upper null hypothesis value for the
#'   info-ratio z-test. Default 5 (handoff-canonical).
#' @param tau_ncp Numeric scalar. Lower null hypothesis value.
#'   Default 2 (handoff-canonical).
#' @param aggregation Character scalar in \code{c("any_ncp",
#'   "majority", "per_k")}. Default \code{"any_ncp"}.
#'
#' @return A list with \code{cp_a} (logical scalar, aggregated),
#'   \code{cp_W} (logical scalar, aggregated), \code{cp_a_per_k}
#'   (logical vector of length \code{p}), \code{cp_W_per_k}
#'   (logical vector of length \code{p}), \code{cs_model_ncp}
#'   (compiled NCP cmdstan model object), and \code{report} (an
#'   object of class \code{gdpar_preflight_report}).
#' @keywords internal
#' @noRd
preflight_parametrization_multi <- function(prior, stan_data, amm,
                                            preflight_seed = NULL,
                                            verbose = FALSE,
                                            tau_cp = 5, tau_ncp = 2,
                                            aggregation = "any_ncp") {
  assert_inherits(prior, "gdpar_prior", "prior")
  assert_inherits(amm, "amm_spec", "amm")
  if (is.null(amm$p) || amm$p < 1L) {
    gdpar_abort(
      "Internal error: amm$p must be set and positive for the multivariate preflight.",
      class = "gdpar_internal_error"
    )
  }
  if (!aggregation %in% c("any_ncp", "majority", "per_k")) {
    gdpar_abort(
      sprintf(
        "Argument 'aggregation' must be one of {'any_ncp', 'majority', 'per_k'}; received '%s'.",
        as.character(aggregation)
      ),
      class = "gdpar_input_error",
      data = list(received = aggregation)
    )
  }

  p <- as.integer(amm$p)
  has_a <- isTRUE(stan_data$use_a == 1L)
  has_W <- isTRUE(stan_data$use_W == 1L)
  J_a_per_k <- as.integer(stan_data$J_a_per_k %||% rep(0L, p))
  J_a_free_per_k <- pmax(0L, J_a_per_k)
  if (!has_a) J_a_free_per_k <- rep(0L, p)
  a_raw_offset_per_k <- c(0L, cumsum(J_a_free_per_k)[-p])
  W_per_k_dim <- as.integer(stan_data$W_per_k_dim %||% 0L)
  d <- as.integer(stan_data$d %||% 0L)
  theta_anchor_vec <- as.double(stan_data$theta_anchor)

  needs_a_per_k <- has_a & (J_a_free_per_k >= 1L)
  needs_W_per_k <- rep(has_W & W_per_k_dim >= 1L & d >= 1L, p)

  per_dim_rows <- vector("list", 2L * p)
  decision_a <- rep("absent", p)
  decision_W <- rep("absent", p)
  reason_a <- rep("absent_or_degenerate", p)
  reason_W <- rep("absent_or_degenerate", p)
  t_attr_a <- rep(NA_real_, p)
  t_attr_W <- rep(NA_real_, p)
  t_info_cp_a <- rep(NA_real_, p)
  t_info_ncp_a <- rep(NA_real_, p)
  t_info_cp_W <- rep(NA_real_, p)
  t_info_ncp_W <- rep(NA_real_, p)
  n_divergent <- NA_integer_
  div_pct <- NA_real_
  ebfmi_min <- NA_real_

  cs_model_ncp <- NULL
  used_preflight <- FALSE

  if (any(needs_a_per_k) || any(needs_W_per_k)) {
    stan_src <- generate_stan_code_multi(prior, cp_a = FALSE, cp_W = FALSE)
    stan_path <- write_stan_to_tempfile(stan_src)
    cs_model_ncp <- cmdstanr::cmdstan_model(stan_path)
    pre_fit <- run_preflight_sample(cs_model_ncp, stan_data,
                                    preflight_seed, verbose)
    used_preflight <- TRUE
    diag <- pre_fit$diagnostic_summary(quiet = TRUE)
    total_transitions <- 2L * 200L
    n_divergent <- as.integer(sum(diag$num_divergent))
    div_pct <- n_divergent / total_transitions
    ebfmi_min <- min(diag$ebfmi)

    alpha <- 0.025
    crit_t <- qt(1 - alpha, df = max(1L, n_divergent - 1L))
    crit_z <- qnorm(1 - alpha)

    # Filter 1 (attribution) per-k --------------------------------
    if (n_divergent >= 1L) {
      for (k in seq_len(p)) {
        if (needs_a_per_k[k]) {
          s <- preflight_attribution_score_a_per_k_from_fit(
            pre_fit, k, J_a_free_per_k[k], a_raw_offset_per_k[k]
          )
          t_attr_a[k] <- s
          if (!is.na(s) && s > crit_t) {
            decision_a[k] <- "CP"
            reason_a[k] <- "filter_attribution"
          }
        }
        if (needs_W_per_k[k]) {
          s <- preflight_attribution_score_W_per_k_from_fit(
            pre_fit, k, W_per_k_dim, d
          )
          t_attr_W[k] <- s
          if (!is.na(s) && s > crit_t) {
            decision_W[k] <- "CP"
            reason_W[k] <- "filter_attribution"
          }
        }
      }
    }

    # Filter 2 (E-BFMI): global, flips undecided to CP ------------
    if (isTRUE(ebfmi_min < 0.3)) {
      for (k in seq_len(p)) {
        if (needs_a_per_k[k] && decision_a[k] == "absent") {
          decision_a[k] <- "CP"
          reason_a[k] <- "filter_ebfmi"
        }
        if (needs_W_per_k[k] && decision_W[k] == "absent") {
          decision_W[k] <- "CP"
          reason_W[k] <- "filter_ebfmi"
        }
      }
    }

    # Filter 3 (info ratio Path B') per-k -------------------------
    for (k in seq_len(p)) {
      if (needs_a_per_k[k] && decision_a[k] == "absent") {
        res <- preflight_a_info_per_k_from_fit(
          pre_fit, k, J_a_free_per_k[k],
          a_raw_offset_per_k[k],
          tau_cp = tau_cp, tau_ncp = tau_ncp
        )
        t_info_cp_a[k] <- res$t_cp
        t_info_ncp_a[k] <- res$t_ncp
        if (is.na(res$t_cp)) {
          decision_a[k] <- "NCP"
          reason_a[k] <- "filter_info_undefined_ncp"
        } else if (res$t_cp > crit_z) {
          decision_a[k] <- "CP"
          reason_a[k] <- "filter_info_high"
        } else if (res$t_ncp < -crit_z) {
          decision_a[k] <- "NCP"
          reason_a[k] <- "filter_info_low"
        } else {
          decision_a[k] <- "NCP"
          reason_a[k] <- "filter_info_ambiguous_ncp"
        }
      } else if (!needs_a_per_k[k]) {
        decision_a[k] <- "absent"
      }
      if (needs_W_per_k[k] && decision_W[k] == "absent") {
        res <- preflight_W_info_per_k_from_fit(
          pre_fit, k, W_per_k_dim, d,
          theta_anchor_k = theta_anchor_vec[k],
          tau_cp = tau_cp, tau_ncp = tau_ncp
        )
        t_info_cp_W[k] <- res$t_cp
        t_info_ncp_W[k] <- res$t_ncp
        if (is.na(res$t_cp)) {
          decision_W[k] <- "NCP"
          reason_W[k] <- "filter_info_undefined_ncp"
        } else if (res$t_cp > crit_z) {
          decision_W[k] <- "CP"
          reason_W[k] <- "filter_info_high"
        } else if (res$t_ncp < -crit_z) {
          decision_W[k] <- "NCP"
          reason_W[k] <- "filter_info_low"
        } else {
          decision_W[k] <- "NCP"
          reason_W[k] <- "filter_info_ambiguous_ncp"
        }
      } else if (!needs_W_per_k[k]) {
        decision_W[k] <- "absent"
      }
    }
  }

  per_dim <- data.frame(
    component       = c(rep("a", p), rep("W", p)),
    dim             = c(seq_len(p), seq_len(p)),
    decision        = c(decision_a, decision_W),
    decision_reason = c(reason_a, reason_W),
    n_divergent     = rep(n_divergent, 2L * p),
    div_pct         = rep(div_pct, 2L * p),
    ebfmi_min       = rep(ebfmi_min, 2L * p),
    t_attribution   = c(t_attr_a, t_attr_W),
    t_info_cp       = c(t_info_cp_a, t_info_cp_W),
    t_info_ncp      = c(t_info_ncp_a, t_info_ncp_W),
    stringsAsFactors = FALSE
  )

  settings <- list(
    tau_cp        = tau_cp,
    tau_ncp       = tau_ncp,
    aggregation   = aggregation,
    used_preflight = used_preflight,
    n_chains      = 2L,
    n_warmup      = 200L,
    n_sampling    = 200L,
    adapt_delta   = 0.95,
    max_treedepth = 10L
  )
  report <- new_gdpar_preflight_report(per_dim, method = aggregation,
                                       settings = settings)
  global <- report$global
  cp_a_global <- decision_to_logical(
    global$global_decision[global$component == "a"]
  )
  cp_W_global <- decision_to_logical(
    global$global_decision[global$component == "W"]
  )
  cp_a_per_k <- decision_to_logical(decision_a)
  cp_W_per_k <- decision_to_logical(decision_W)
  list(
    cp_a         = cp_a_global,
    cp_W         = cp_W_global,
    cp_a_per_k   = cp_a_per_k,
    cp_W_per_k   = cp_W_per_k,
    cs_model_ncp = cs_model_ncp,
    report       = report
  )
}


#' Internal: per-k Path B' info-ratio extractor for component a
#'
#' Extracts the per-draw chunk of \code{a_coef[k][1:J_a_free[k]]}
#' (the free additive coefficients) and the per-draw \code{sigma_a[k]}
#' from the preflight fit and delegates to \code{preflight_info_ratio_t}.
#'
#' @keywords internal
#' @noRd
preflight_a_info_per_k_from_fit <- function(fit, k, J_a_free_k,
                                            a_raw_offset_k,
                                            tau_cp, tau_ncp) {
  if (J_a_free_k < 1L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  draws_sigma <- tryCatch(
    fit$draws(variables = sprintf("sigma_a[%d]", k),
              format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_sigma) || ncol(draws_sigma) == 0L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  sigma_vec <- as.numeric(draws_sigma[, 1])
  if (any(is.na(sigma_vec)) || any(sigma_vec <= 0)) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  vars <- sprintf("a_coef[%d,%d]", k, seq_len(J_a_free_k))
  draws_coef <- tryCatch(
    fit$draws(variables = vars, format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_coef) || ncol(draws_coef) != J_a_free_k) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  coef_mat <- as.matrix(draws_coef)
  preflight_info_ratio_t(
    effective_coef = coef_mat,
    reference_scale_per_draw = sigma_vec,
    n_chains = 2L,
    tau_cp = tau_cp, tau_ncp = tau_ncp
  )
}


#' Internal: per-k Path B' info-ratio extractor for component W
#'
#' Reconstructs the per-draw effective coefficient at coordinate
#' \eqn{jj} for the k-th coord of theta_ref as
#' \deqn{ \mathrm{eff}[t, jj] = \sum_{m=1}^{W_{per_k_dim}}
#' (\theta_{ref}[t, k]^m - \theta_{anchor}[k]^m)
#' W_{raw}[t, (k-1) W_{per_k_dim} + m, jj] \sigma_W[t] } and the
#' per-draw reference scale as \eqn{\sigma_W[t] \sqrt{\sum_m
#' (\theta_{ref}[t, k]^m - \theta_{anchor}[k]^m)^2}}, then delegates
#' to \code{preflight_info_ratio_t}.
#'
#' @keywords internal
#' @noRd
preflight_W_info_per_k_from_fit <- function(fit, k, W_per_k_dim, d,
                                            theta_anchor_k,
                                            tau_cp, tau_ncp) {
  if (W_per_k_dim < 1L || d < 1L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  draws_theta <- tryCatch(
    fit$draws(variables = sprintf("theta_ref[%d]", k),
              format = "draws_matrix"),
    error = function(e) NULL
  )
  draws_sigma_W <- tryCatch(
    fit$draws(variables = "sigma_W[1]", format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_theta) || ncol(draws_theta) == 0L ||
      is.null(draws_sigma_W) || ncol(draws_sigma_W) == 0L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  theta_vec <- as.numeric(draws_theta[, 1])
  sigma_vec <- as.numeric(draws_sigma_W[, 1])
  if (any(is.na(theta_vec)) || any(is.na(sigma_vec)) ||
      any(sigma_vec <= 0)) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  n_draws <- length(theta_vec)
  block_rows <- (k - 1L) * W_per_k_dim + seq_len(W_per_k_dim)
  vars <- as.vector(outer(block_rows, seq_len(d),
                          function(r, c) sprintf("W_raw[%d,%d]", r, c)))
  draws_W <- tryCatch(
    fit$draws(variables = vars, format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_W) || ncol(draws_W) != W_per_k_dim * d) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  W_block_mat <- as.matrix(draws_W)
  W_block_arr <- array(W_block_mat, dim = c(n_draws, W_per_k_dim, d))
  basis_diff <- vapply(
    seq_len(W_per_k_dim),
    function(m) theta_vec^m - theta_anchor_k^m,
    numeric(n_draws)
  )
  if (!is.matrix(basis_diff)) {
    basis_diff <- matrix(basis_diff, nrow = n_draws,
                         ncol = W_per_k_dim)
  }
  basis_norm <- sqrt(rowSums(basis_diff^2))
  if (any(!is.finite(basis_norm)) || any(basis_norm <= 0)) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  reference_scale_per_draw <- sigma_vec * basis_norm
  eff_coef <- matrix(0, nrow = n_draws, ncol = d)
  for (jj in seq_len(d)) {
    contrib <- numeric(n_draws)
    for (m in seq_len(W_per_k_dim)) {
      contrib <- contrib + basis_diff[, m] *
        W_block_arr[, m, jj] * sigma_vec
    }
    eff_coef[, jj] <- contrib
  }
  preflight_info_ratio_t(
    effective_coef = eff_coef,
    reference_scale_per_draw = reference_scale_per_draw,
    n_chains = 2L,
    tau_cp = tau_cp, tau_ncp = tau_ncp
  )
}


#' Internal: per-k attribution score for component a (Path B' filter 1)
#'
#' Extracts \code{sigma_a[k]} and the flat-packed chunk
#' \code{a_raw[a_raw_offset_k + 1:J_a_free_k]} from the preflight fit
#' and delegates to \code{preflight_attribution_score}.
#'
#' @keywords internal
#' @noRd
preflight_attribution_score_a_per_k_from_fit <- function(fit, k,
                                                         J_a_free_k,
                                                         a_raw_offset_k) {
  if (J_a_free_k < 1L) return(NA_real_)
  draws_sigma <- tryCatch(
    fit$draws(variables = sprintf("sigma_a[%d]", k),
              format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_sigma) || ncol(draws_sigma) == 0L) return(NA_real_)
  sigma_vec <- as.numeric(draws_sigma[, 1])
  if (any(sigma_vec <= 0) || any(is.na(sigma_vec))) return(NA_real_)
  vars <- sprintf("a_raw[%d]",
                  a_raw_offset_k + seq_len(J_a_free_k))
  draws_raw <- tryCatch(
    fit$draws(variables = vars, format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_raw) || ncol(draws_raw) != J_a_free_k) {
    return(NA_real_)
  }
  raw_mat <- as.matrix(draws_raw)
  raw_norm <- sqrt(rowSums(raw_mat^2))
  diag_mat <- tryCatch(
    fit$sampler_diagnostics(format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(diag_mat) || !"divergent__" %in% colnames(diag_mat)) {
    return(NA_real_)
  }
  divergent <- as.logical(diag_mat[, "divergent__"])
  preflight_attribution_score(log(sigma_vec), raw_norm, divergent)
}


#' Internal: per-k attribution score for component W (Path B' filter 1)
#'
#' Extracts \code{sigma_W[1]} (global) and the rows of \code{W_raw}
#' that belong to the k-th block, computes the per-draw norm of that
#' block, and delegates to \code{preflight_attribution_score}.
#'
#' @keywords internal
#' @noRd
preflight_attribution_score_W_per_k_from_fit <- function(fit, k,
                                                         W_per_k_dim,
                                                         d) {
  if (W_per_k_dim < 1L || d < 1L) return(NA_real_)
  draws_sigma <- tryCatch(
    fit$draws(variables = "sigma_W[1]", format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_sigma) || ncol(draws_sigma) == 0L) return(NA_real_)
  sigma_vec <- as.numeric(draws_sigma[, 1])
  if (any(sigma_vec <= 0) || any(is.na(sigma_vec))) return(NA_real_)
  block_rows <- (k - 1L) * W_per_k_dim + seq_len(W_per_k_dim)
  vars <- as.vector(outer(block_rows, seq_len(d),
                          function(r, c) sprintf("W_raw[%d,%d]", r, c)))
  draws_W <- tryCatch(
    fit$draws(variables = vars, format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_W) || ncol(draws_W) != W_per_k_dim * d) {
    return(NA_real_)
  }
  W_block_mat <- as.matrix(draws_W)
  W_norm <- sqrt(rowSums(W_block_mat^2))
  diag_mat <- tryCatch(
    fit$sampler_diagnostics(format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(diag_mat) || !"divergent__" %in% colnames(diag_mat)) {
    return(NA_real_)
  }
  divergent <- as.logical(diag_mat[, "divergent__"])
  preflight_attribution_score(log(sigma_vec), W_norm, divergent)
}


#' Internal: map decision string(s) to a logical CP flag
#'
#' Each element of \code{x} is mapped as follows: "CP" -> TRUE,
#' "NCP" -> FALSE, "absent" -> FALSE, "per_k" -> NA (sentinel; the
#' resolver translates this to an unsupported-feature error in
#' Phase H.1).
#'
#' @keywords internal
#' @noRd
decision_to_logical <- function(x) {
  vapply(x, function(v) {
    switch(v,
      "CP"     = TRUE,
      "NCP"    = FALSE,
      "absent" = FALSE,
      "per_k"  = NA,
      gdpar_abort(
        sprintf("Unknown decision value '%s'.", v),
        class = "gdpar_internal_error",
        data = list(value = v)
      )
    )
  }, logical(1L), USE.NAMES = FALSE)
}


#' Internal: NULL-coalescing operator (local)
#'
#' Local fallback so the helpers above do not depend on rlang or
#' base-R \code{\%||\%} availability across R versions.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Internal: resolve parametrization flags for the long fit (multi)
#'
#' Multivariate counterpart of \code{resolve_parametrization()}.
#' Combines the user-supplied \code{parametrization} (global),
#' \code{parametrization_a} / \code{parametrization_W}
#' (per-component, scalar), and the optional pre-flight diagnostic.
#' The per-k toggle is not exposed to user input in Phase H.1;
#' \code{aggregation = "per_k"} is rejected with
#' \code{gdpar_unsupported_feature_error} because the Stan template
#' wiring for per-k priors is delivered in Phase H.2.
#'
#' @inheritParams preflight_parametrization_multi
#' @param parametrization Character: \code{"auto"}, \code{"ncp"} or
#'   \code{"cp"} (already validated by the caller).
#' @param parametrization_a \code{NULL} or \code{"ncp"}. \code{"cp"}
#'   is rejected with \code{gdpar_unsupported_feature_error} until
#'   Phase H.2 lands.
#' @param parametrization_W \code{NULL} or \code{"ncp"}/\code{"cp"}.
#'   \code{"cp"} is honored (it sets \code{cp_W = TRUE} for every
#'   coordinate via the global \code{sigma_W[1]}).
#'
#' @return List with \code{cp_a} (logical), \code{cp_W} (logical),
#'   \code{cp_a_per_k} (logical vector length p),
#'   \code{cp_W_per_k} (logical vector length p), \code{report}
#'   (\code{gdpar_preflight_report} or \code{NULL} when the
#'   pre-flight is skipped).
#' @keywords internal
#' @noRd
resolve_parametrization_multi <- function(parametrization,
                                          parametrization_a,
                                          parametrization_W,
                                          parametrization_aggregation,
                                          prior, stan_data, amm,
                                          preflight_seed,
                                          verbose,
                                          tau_cp = 5, tau_ncp = 2) {
  p <- as.integer(amm$p)
  if (is.null(parametrization_aggregation)) {
    parametrization_aggregation <- "any_ncp"
  }
  if (!parametrization_aggregation %in%
      c("any_ncp", "majority", "per_k")) {
    gdpar_abort(
      sprintf(
        "Argument 'parametrization_aggregation' must be one of {'any_ncp', 'majority', 'per_k'}; received '%s'.",
        as.character(parametrization_aggregation)
      ),
      class = "gdpar_input_error",
      data = list(received = parametrization_aggregation)
    )
  }

  cp_a_user <- if (!is.null(parametrization_a)) {
    parametrization_a == "cp"
  } else NA
  cp_W_user <- if (!is.null(parametrization_W)) {
    parametrization_W == "cp"
  } else NA
  cp_a_global <- switch(parametrization,
                        ncp = FALSE, cp = TRUE, auto = NA)
  cp_W_global <- cp_a_global
  needs_pf_a <- is.na(cp_a_user) && is.na(cp_a_global)
  needs_pf_W <- is.na(cp_W_user) && is.na(cp_W_global)

  if (needs_pf_a || needs_pf_W) {
    pf <- preflight_parametrization_multi(
      prior, stan_data, amm,
      preflight_seed = preflight_seed,
      verbose = verbose,
      tau_cp = tau_cp, tau_ncp = tau_ncp,
      aggregation = parametrization_aggregation
    )
    cp_a <- if (!is.na(cp_a_user)) cp_a_user else pf$cp_a
    cp_W <- if (!is.na(cp_W_user)) cp_W_user else pf$cp_W
    cp_a_per_k <- if (!is.na(cp_a_user)) {
      rep(cp_a_user, p)
    } else pf$cp_a_per_k
    cp_W_per_k <- if (!is.na(cp_W_user)) {
      rep(cp_W_user, p)
    } else pf$cp_W_per_k
    report <- pf$report
  } else {
    cp_a <- if (!is.na(cp_a_user)) cp_a_user else cp_a_global
    cp_W <- if (!is.na(cp_W_user)) cp_W_user else cp_W_global
    cp_a_per_k <- rep(cp_a, p)
    cp_W_per_k <- rep(cp_W, p)
    report <- NULL
  }

  # Phase H.2 contract: the multivariate Stan template honors per-k
  # CP/NCP for the additive component 'a' via segment()-based priors
  # (uniform branches preserve bit-exact paridad with H.1 NCP-only).
  # The cp_a global is recomputed from cp_a_per_k under the chosen
  # aggregation strategy when "per_k" is requested it stays a vector;
  # otherwise it collapses to a scalar that drives the uniform
  # branch.
  cp_a_resolved <- cp_a
  cp_a_per_k_resolved <- cp_a_per_k
  if (identical(parametrization_aggregation, "per_k")) {
    cp_a_resolved <- NA
  } else {
    if (length(unique(cp_a_per_k_resolved)) > 1L) {
      cp_a_per_k_resolved <- rep(cp_a_resolved, p)
    }
  }

  # Mirror the resolution for cp_W. The current model uses a single
  # global sigma_W[1] shared across blocks, so per-k W decisions can
  # be honored at the sampling level only by aggregating to a global
  # cp_W. When per_k aggregation is requested but the preflight per-k
  # for W is heterogeneous, emit a structured informational message
  # documenting the limitation; the report retains the per-k detail.
  cp_W_resolved <- cp_W
  cp_W_per_k_resolved <- cp_W_per_k
  if (length(unique(cp_W_per_k_resolved)) > 1L) {
    gdpar_inform(
      paste(
        "preflight detected heterogeneous per-k decisions for W, but",
        "the current model uses a single global sigma_W[1] shared",
        "across blocks. The W per-k decisions are recorded in",
        "fit$parametrization$report for auditability; the sampler",
        "honors only the aggregated cp_W. Bloque 8 (multi-parametric",
        "extension) may promote sigma_W to array[p] and honor per-k W."
      ),
      class = "gdpar_W_per_k_heterogeneous_message",
      data = list(cp_W_per_k = cp_W_per_k_resolved)
    )
    cp_W_per_k_resolved <- rep(cp_W_resolved, p)
  }

  list(
    cp_a         = as.logical(cp_a_resolved),
    cp_W         = as.logical(cp_W_resolved),
    cp_a_per_k   = as.logical(cp_a_per_k_resolved),
    cp_W_per_k   = as.logical(cp_W_per_k_resolved),
    aggregation  = parametrization_aggregation,
    report       = report
  )
}

#' Internal: pre-flight diagnostic to select CP or NCP per component
#'
#' Compiles and samples a short non-centered (NCP) version of the
#' Path 1 Stan model and applies three filters to decide whether the
#' additive component \code{a} and/or the modulating component
#' \code{W} should switch to the centered parametrization (CP) for
#' the long fit. The decision is dataset-driven: filters 4 and 5
#' use t-statistics against closed-form null distributions; filter
#' 6 uses an asymptotic z-statistic with a chain-aware
#' block-bootstrap standard error (Path B', see below).
#'
#' Filters, in priority order:
#' \enumerate{
#'   \item Divergence attribution (only if \eqn{n_{div} \ge 1}). For
#'         each component, score \eqn{S = \mathrm{mean}_{div}(-z_{\log
#'         \sigma} z_{||raw||})} where the z-scores are computed over
#'         all pre-flight transitions. Center by the global posterior
#'         mean of the score so \eqn{H_0: E[S_{centered}] = 0}.
#'         Standard error \eqn{SE = \sqrt{\mathrm{Var}(g) / n_{div}}}
#'         with \eqn{\mathrm{Var}(g)} from all transitions. Reject at
#'         \eqn{\alpha = 0.025} one-sided.
#'   \item E-BFMI: if \code{ebfmi_min < 0.3}, force CP on any
#'         component still undecided. Energy mixing pathology is not
#'         attributable per-component.
#'   \item Info ratio (Path B'). Build a per-draw effective
#'         coefficient \eqn{\mathrm{eff}[t, jj]} and a per-draw
#'         reference scale \eqn{\mathrm{ref}[t]} such that
#'         \eqn{\mathrm{ref}[t]} is the conditional prior standard
#'         deviation of \eqn{\mathrm{eff}[t, jj]} given the
#'         hyperparameters at draw \eqn{t}. Compute the per-coordinate
#'         log info ratio \eqn{\log(\mathrm{mean}(\mathrm{ref}) /
#'         \mathrm{sd}(\mathrm{eff}_j))} on the full data and average
#'         across coordinates; estimate the standard error of that
#'         mean by a chain-aware block bootstrap (block size 10) to
#'         mitigate MCMC autocorrelation. Two asymptotic z-tests
#'         against null hypotheses \eqn{\log \tau_{cp}} (upper) and
#'         \eqn{\log \tau_{ncp}} (lower). Reject toward CP if the
#'         upper test passes; toward NCP if the lower test passes;
#'         otherwise NCP conservative. The effective-weight
#'         construction handles the rank deficiency of the
#'         per-coordinate posterior in \eqn{W} when the polynomial
#'         basis is evaluated at a tightly identified scalar
#'         \eqn{\theta_{ref}}; the block bootstrap handles the
#'         within-chain autocorrelation of HMC draws.
#' }
#'
#' Pre-flight settings are fixed (\code{adapt_delta = 0.95},
#' \code{max_treedepth = 10}, 200 warmup + 200 sampling, 2 chains)
#' so the decision depends on the model and data, not on the user's
#' settings for the long fit.
#'
#' @param prior A \code{gdpar_prior} object.
#' @param stan_data List assembled by \code{assemble_stan_data}.
#' @param amm An \code{amm_spec} object.
#' @param preflight_seed Integer or \code{NULL}; forwarded to the
#'   pre-flight sampler.
#' @param verbose Logical. If \code{FALSE}, pre-flight stderr is
#'   sunk to a tempfile and discarded.
#' @param tau_cp Numeric scalar. Upper null hypothesis value for the
#'   info-ratio z-test (filter 6, Path B'). Defaults to 5,
#'   calibrated empirically against the eight canonical scenarios in
#'   \code{inst/benchmarks/calibrate_cp_ncp.R} (overall hit rate
#'   0.92, individual hit rate 1.00 on \code{a} and 0.83 on
#'   \code{W}, with the single residual miss on \code{W}
#'   attributable to structural confounding when \code{a} and
#'   \code{W} share covariates).
#' @param tau_ncp Numeric scalar. Lower null hypothesis value for
#'   the info-ratio z-test. Defaults to 2; same calibration.
#' @return A list with elements \code{cp_a} (logical), \code{cp_W}
#'   (logical), \code{cs_model_ncp} (compiled NCP cmdstan model
#'   object, suitable for reuse if both decisions are NCP), and
#'   \code{meta} (diagnostic statistics and the textual reason per
#'   component).
#' @keywords internal
#' @noRd
preflight_parametrization <- function(prior, stan_data, amm,
                                      preflight_seed = NULL,
                                      verbose = FALSE,
                                      tau_cp = 5, tau_ncp = 2) {
  assert_inherits(prior, "gdpar_prior", "prior")
  assert_inherits(amm, "amm_spec", "amm")

  has_a <- !is.null(amm$a)
  has_W <- !is.null(amm$W)
  J_a <- as.integer(stan_data$J_a)
  J_a_free <- if (has_a && J_a >= 1L) J_a else 0L
  dim_W <- as.integer(stan_data$dim_W)
  d <- as.integer(stan_data$d)
  needs_a <- has_a && J_a_free >= 1L
  needs_W <- has_W && dim_W >= 1L && d >= 1L
  n_coords_W <- as.integer(dim_W * d)

  meta <- list(
    used_preflight       = FALSE,
    n_divergent          = NA_integer_,
    div_pct              = NA_real_,
    ebfmi_min            = NA_real_,
    t_attribution_a      = NA_real_,
    t_attribution_W      = NA_real_,
    t_info_cp_a          = NA_real_,
    t_info_ncp_a         = NA_real_,
    t_info_cp_W          = NA_real_,
    t_info_ncp_W         = NA_real_,
    decision_reason_a    = if (needs_a) "" else "absent_or_degenerate",
    decision_reason_W    = if (needs_W) "" else "absent_or_degenerate"
  )

  if (!needs_a && !needs_W) {
    return(list(
      cp_a = FALSE, cp_W = FALSE,
      cs_model_ncp = NULL, meta = meta
    ))
  }

  stan_src <- generate_stan_code(prior, cp_a = FALSE, cp_W = FALSE)
  stan_path <- write_stan_to_tempfile(stan_src)
  cs_model_ncp <- cmdstanr::cmdstan_model(stan_path)

  pre_fit <- run_preflight_sample(cs_model_ncp, stan_data,
                                  preflight_seed, verbose)
  meta$used_preflight <- TRUE

  diag <- pre_fit$diagnostic_summary(quiet = TRUE)
  total_transitions <- 2L * 200L
  n_div <- as.integer(sum(diag$num_divergent))
  meta$n_divergent <- n_div
  meta$div_pct <- n_div / total_transitions
  meta$ebfmi_min <- min(diag$ebfmi)

  alpha <- 0.025
  cp_a <- NA
  cp_W <- NA

  if (needs_a && n_div >= 1L) {
    score_a <- preflight_attribution_score_from_fit(pre_fit, "a")
    meta$t_attribution_a <- score_a
    crit <- qt(1 - alpha, df = max(1L, n_div - 1L))
    if (!is.na(score_a) && score_a > crit) {
      cp_a <- TRUE
      meta$decision_reason_a <- "filter_attribution"
    }
  }
  if (needs_W && n_div >= 1L) {
    score_W <- preflight_attribution_score_from_fit(pre_fit, "W")
    meta$t_attribution_W <- score_W
    crit <- qt(1 - alpha, df = max(1L, n_div - 1L))
    if (!is.na(score_W) && score_W > crit) {
      cp_W <- TRUE
      meta$decision_reason_W <- "filter_attribution"
    }
  }

  if (meta$ebfmi_min < 0.3) {
    if (needs_a && is.na(cp_a)) {
      cp_a <- TRUE
      meta$decision_reason_a <- "filter_ebfmi"
    }
    if (needs_W && is.na(cp_W)) {
      cp_W <- TRUE
      meta$decision_reason_W <- "filter_ebfmi"
    }
  }

  crit_z <- qnorm(1 - alpha)

  if (needs_a && is.na(cp_a)) {
    res_a <- preflight_info_ratio_t_from_fit(
      pre_fit, "a", J_a_free,
      n_chains = 2L,
      tau_cp = tau_cp, tau_ncp = tau_ncp
    )
    meta$t_info_cp_a <- res_a$t_cp
    meta$t_info_ncp_a <- res_a$t_ncp
    if (is.na(res_a$t_cp)) {
      cp_a <- FALSE
      meta$decision_reason_a <- "filter_info_undefined_ncp"
    } else if (res_a$t_cp > crit_z) {
      cp_a <- TRUE
      meta$decision_reason_a <- "filter_info_high"
    } else if (res_a$t_ncp < -crit_z) {
      cp_a <- FALSE
      meta$decision_reason_a <- "filter_info_low"
    } else {
      cp_a <- FALSE
      meta$decision_reason_a <- "filter_info_ambiguous_ncp"
    }
  }

  if (needs_W && is.na(cp_W)) {
    res_W <- preflight_info_ratio_t_from_fit(
      pre_fit, "W", n_coords_W,
      n_chains = 2L,
      theta_anchor = stan_data$theta_anchor,
      dim_W = dim_W, d = d,
      tau_cp = tau_cp, tau_ncp = tau_ncp
    )
    meta$t_info_cp_W <- res_W$t_cp
    meta$t_info_ncp_W <- res_W$t_ncp
    if (is.na(res_W$t_cp)) {
      cp_W <- FALSE
      meta$decision_reason_W <- "filter_info_undefined_ncp"
    } else if (res_W$t_cp > crit_z) {
      cp_W <- TRUE
      meta$decision_reason_W <- "filter_info_high"
    } else if (res_W$t_ncp < -crit_z) {
      cp_W <- FALSE
      meta$decision_reason_W <- "filter_info_low"
    } else {
      cp_W <- FALSE
      meta$decision_reason_W <- "filter_info_ambiguous_ncp"
    }
  }

  if (is.na(cp_a)) cp_a <- FALSE
  if (is.na(cp_W)) cp_W <- FALSE

  list(
    cp_a = as.logical(cp_a),
    cp_W = as.logical(cp_W),
    cs_model_ncp = cs_model_ncp,
    meta = meta
  )
}

#' Internal: run the pre-flight sample with stderr handling
#'
#' Wraps \code{cs_model$sample()} with the fixed pre-flight settings
#' (200 warmup + 200 sampling, 2 chains, \code{adapt_delta = 0.95},
#' \code{max_treedepth = 10}). When \code{verbose = FALSE}, stderr is
#' redirected to a tempfile that is deleted on exit so the caller's
#' console stays clean. Restoration of the sink is guaranteed by
#' \code{on.exit} even if the sample throws.
#'
#' @keywords internal
#' @noRd
run_preflight_sample <- function(cs_model, stan_data, seed, verbose) {
  args <- list(
    data            = stan_data,
    chains          = 2L,
    iter_warmup     = 200L,
    iter_sampling   = 200L,
    adapt_delta     = 0.95,
    max_treedepth   = 10L,
    refresh         = 0L,
    show_messages   = verbose,
    show_exceptions = verbose
  )
  if (!is.null(seed)) args$seed <- as.integer(seed)
  if (isTRUE(verbose)) {
    do.call(cs_model$sample, args)
  } else {
    msg_file <- tempfile(fileext = ".log")
    msg_con <- file(msg_file, open = "wt")
    sink(msg_con, type = "message")
    on.exit({
      if (sink.number(type = "message") > 0L) {
        sink(NULL, type = "message")
      }
      if (isOpen(msg_con)) close(msg_con)
      unlink(msg_file)
    }, add = TRUE, after = FALSE)
    do.call(cs_model$sample, args)
  }
}

#' Internal: extract attribution inputs from the pre-flight fit
#'
#' Pulls draws of \code{sigma_X[1]}, \code{X_raw}, and the divergent
#' flag for the requested component, then delegates the t-statistic
#' computation to \code{preflight_attribution_score()}.
#'
#' @keywords internal
#' @noRd
preflight_attribution_score_from_fit <- function(fit, component) {
  if (component == "a") {
    sigma_var <- "sigma_a"
    raw_var <- "a_raw"
  } else if (component == "W") {
    sigma_var <- "sigma_W"
    raw_var <- "W_raw"
  } else {
    return(NA_real_)
  }
  draws_sigma <- tryCatch(
    fit$draws(variables = sigma_var, format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_sigma) || ncol(draws_sigma) == 0L) {
    return(NA_real_)
  }
  draws_raw <- tryCatch(
    fit$draws(variables = raw_var, format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(draws_raw) || ncol(draws_raw) == 0L) {
    return(NA_real_)
  }
  diag_mat <- tryCatch(
    fit$sampler_diagnostics(format = "draws_matrix"),
    error = function(e) NULL
  )
  if (is.null(diag_mat) || !"divergent__" %in% colnames(diag_mat)) {
    return(NA_real_)
  }
  sigma_vec <- as.numeric(draws_sigma[, 1])
  if (any(sigma_vec <= 0) || any(is.na(sigma_vec))) {
    return(NA_real_)
  }
  raw_mat <- as.matrix(draws_raw)
  raw_norm <- sqrt(rowSums(raw_mat^2))
  divergent <- as.logical(diag_mat[, "divergent__"])
  preflight_attribution_score(log(sigma_vec), raw_norm, divergent)
}

#' Internal: pure t-statistic for divergence attribution
#'
#' Inputs are aligned vectors over the flat sequence of pre-flight
#' transitions. The score per transition is \eqn{g_t = -z_{\log
#' \sigma}(t) \cdot z_{||raw||}(t)}; high values flag the funnel
#' region. The test statistic is the centered mean over divergent
#' transitions divided by \eqn{SE = \sqrt{\mathrm{Var}(g) / n_{div}}}
#' with \eqn{\mathrm{Var}(g)} from all transitions (under \eqn{H_0}
#' the divergent subsample is random).
#'
#' Returns \code{NA_real_} when the test cannot be computed (zero
#' variance, no divergences, or length mismatch).
#'
#' @keywords internal
#' @noRd
preflight_attribution_score <- function(log_sigma, raw_norm, divergent) {
  n <- length(log_sigma)
  if (n != length(raw_norm) || n != length(divergent)) {
    return(NA_real_)
  }
  sd_ls <- sd(log_sigma)
  sd_rn <- sd(raw_norm)
  if (is.na(sd_ls) || sd_ls <= 0 || is.na(sd_rn) || sd_rn <= 0) {
    return(NA_real_)
  }
  z_ls <- (log_sigma - mean(log_sigma)) / sd_ls
  z_rn <- (raw_norm - mean(raw_norm)) / sd_rn
  g <- -z_ls * z_rn
  g_div <- g[divergent]
  n_div <- length(g_div)
  if (n_div < 1L) return(NA_real_)
  s_div <- mean(g_div)
  s_centered <- s_div - mean(g)
  sd_g <- sd(g)
  if (is.na(sd_g) || sd_g <= 0) return(NA_real_)
  se <- sd_g / sqrt(n_div)
  s_centered / se
}

#' Internal: chain-aware block bootstrap indices
#'
#' Returns an integer vector of length \code{n_draws} of indices
#' into the original draws matrix, formed by resampling contiguous
#' blocks of \code{block_size} draws within each chain (no
#' resampling across chains). Mitigates MCMC autocorrelation that
#' would lead an iid bootstrap to underestimate the standard error
#' of any draw-level statistic.
#'
#' Assumes the rows of the draws matrix are ordered chain-first
#' (chain 1 iter 1, ..., chain 1 iter \code{n_draws / n_chains},
#' chain 2 iter 1, ...).
#'
#' @keywords internal
#' @noRd
block_bootstrap_indices <- function(n_draws, n_chains,
                                    block_size = 10L) {
  n_draws <- as.integer(n_draws)
  n_chains <- as.integer(n_chains)
  block_size <- as.integer(block_size)
  if (n_draws %% n_chains != 0L) {
    gdpar_abort(
      "n_draws must be divisible by n_chains.",
      class = "gdpar_input_error",
      data = list(n_draws = n_draws, n_chains = n_chains)
    )
  }
  draws_per_chain <- n_draws %/% n_chains
  if (block_size < 1L || block_size > draws_per_chain) {
    gdpar_abort(
      "block_size must be in [1, n_draws / n_chains].",
      class = "gdpar_input_error",
      data = list(block_size = block_size, draws_per_chain = draws_per_chain)
    )
  }
  blocks_per_chain <- as.integer(ceiling(draws_per_chain / block_size))
  out <- integer(n_chains * draws_per_chain)
  for (cc in seq_len(n_chains)) {
    chain_start <- (cc - 1L) * draws_per_chain + 1L
    n_starts <- draws_per_chain - block_size + 1L
    starts <- sample.int(n_starts, blocks_per_chain, replace = TRUE)
    chain_idx <- integer(blocks_per_chain * block_size)
    for (b in seq_len(blocks_per_chain)) {
      pos <- (b - 1L) * block_size + 1L
      chain_idx[pos:(pos + block_size - 1L)] <-
        chain_start + starts[b] - 1L + seq_len(block_size) - 1L
    }
    chain_idx <- chain_idx[seq_len(draws_per_chain)]
    out_offset <- (cc - 1L) * draws_per_chain
    out[(out_offset + 1L):(out_offset + draws_per_chain)] <- chain_idx
  }
  out
}

#' Internal: extract effective coefficients and reference scale
#' from the pre-flight fit (Path B')
#'
#' For component \code{a}, the effective coefficient is
#' \code{a_coef[, 1:J_a_free]} and the reference scale per draw is
#' \code{sigma_a[1]}. For component \code{W}, the effective
#' coefficient at coordinate \eqn{jj} per draw \eqn{t} is
#' \eqn{\mathrm{eff}[t, jj] = \sum_k (\theta_{ref}[t]^k -
#' \theta_{anchor}^k) W_{raw}[t, k, jj] \sigma_W[t]}, and the
#' reference scale per draw is \eqn{\sigma_W[t] \sqrt{\sum_k
#' (\theta_{ref}[t]^k - \theta_{anchor}^k)^2}} (the conditional
#' prior standard deviation of \eqn{\mathrm{eff}[t, jj]} given the
#' hyperparameters at draw \eqn{t}).
#'
#' @param fit A cmdstanr fit object from the NCP pre-flight.
#' @param component Either \code{"a"} or \code{"W"}.
#' @param n_coords Integer. Used only for component \code{"a"}
#'   (\code{J_a_free}). Ignored for \code{"W"}.
#' @param n_chains Integer.
#' @param theta_anchor Numeric scalar. Required when
#'   \code{component == "W"}.
#' @param dim_W,d Integer. Required when \code{component == "W"}.
#' @param tau_cp,tau_ncp Numeric scalars. Forwarded.
#' @return List with \code{t_cp} and \code{t_ncp}.
#' @keywords internal
#' @noRd
preflight_info_ratio_t_from_fit <- function(fit, component, n_coords,
                                            n_chains,
                                            theta_anchor = NULL,
                                            dim_W = NULL, d = NULL,
                                            tau_cp = 5,
                                            tau_ncp = 2) {
  if (component == "a") {
    draws_coef <- tryCatch(
      fit$draws(variables = "a_coef", format = "draws_matrix"),
      error = function(e) NULL
    )
    draws_sigma <- tryCatch(
      fit$draws(variables = "sigma_a", format = "draws_matrix"),
      error = function(e) NULL
    )
    if (is.null(draws_coef) || ncol(draws_coef) == 0L ||
        is.null(draws_sigma) || ncol(draws_sigma) == 0L) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    coef_mat <- as.matrix(draws_coef)[, seq_len(n_coords),
                                      drop = FALSE]
    sigma_vec <- as.numeric(draws_sigma[, 1])
    if (any(is.na(sigma_vec)) || any(sigma_vec <= 0)) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    return(preflight_info_ratio_t(
      effective_coef = coef_mat,
      reference_scale_per_draw = sigma_vec,
      n_chains = n_chains,
      tau_cp = tau_cp, tau_ncp = tau_ncp
    ))
  }
  if (component == "W") {
    if (is.null(theta_anchor) || is.null(dim_W) || is.null(d)) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    dim_W <- as.integer(dim_W)
    d <- as.integer(d)
    draws_theta <- tryCatch(
      fit$draws(variables = "theta_ref", format = "draws_matrix"),
      error = function(e) NULL
    )
    draws_sigma <- tryCatch(
      fit$draws(variables = "sigma_W", format = "draws_matrix"),
      error = function(e) NULL
    )
    draws_W <- tryCatch(
      fit$draws(variables = "W_raw", format = "draws_matrix"),
      error = function(e) NULL
    )
    if (is.null(draws_theta) || ncol(draws_theta) == 0L ||
        is.null(draws_sigma) || ncol(draws_sigma) == 0L ||
        is.null(draws_W) || ncol(draws_W) == 0L) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    theta_vec <- as.numeric(draws_theta[, 1])
    sigma_vec <- as.numeric(draws_sigma[, 1])
    if (any(is.na(theta_vec)) || any(is.na(sigma_vec)) ||
        any(sigma_vec <= 0)) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    n_draws <- length(theta_vec)
    if (ncol(draws_W) != dim_W * d) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    W_raw_mat <- as.matrix(draws_W)
    W_raw_array <- array(W_raw_mat, dim = c(n_draws, dim_W, d))
    basis_diff <- vapply(seq_len(dim_W),
                         function(k) theta_vec^k - theta_anchor^k,
                         numeric(n_draws))
    if (!is.matrix(basis_diff)) {
      basis_diff <- matrix(basis_diff, nrow = n_draws, ncol = dim_W)
    }
    basis_norm <- sqrt(rowSums(basis_diff^2))
    if (any(!is.finite(basis_norm)) || any(basis_norm <= 0)) {
      return(list(t_cp = NA_real_, t_ncp = NA_real_))
    }
    reference_scale_per_draw <- sigma_vec * basis_norm
    eff_coef <- matrix(0, nrow = n_draws, ncol = d)
    for (jj in seq_len(d)) {
      contrib <- numeric(n_draws)
      for (k in seq_len(dim_W)) {
        contrib <- contrib + basis_diff[, k] *
          W_raw_array[, k, jj] * sigma_vec
      }
      eff_coef[, jj] <- contrib
    }
    return(preflight_info_ratio_t(
      effective_coef = eff_coef,
      reference_scale_per_draw = reference_scale_per_draw,
      n_chains = n_chains,
      tau_cp = tau_cp, tau_ncp = tau_ncp
    ))
  }
  list(t_cp = NA_real_, t_ncp = NA_real_)
}

#' Internal: pure z-statistics for the info-ratio filter (Path B')
#'
#' Computes per coordinate \eqn{\log(\mathrm{mean}(\mathrm{ref}) /
#' \mathrm{sd}(\mathrm{eff}_j))} on the full data and averages
#' across coordinates; estimates the standard error of that mean by
#' a chain-aware block bootstrap of the rows (block size
#' \code{block_size}, \code{n_boot} replicates). Two asymptotic
#' z-tests against null hypotheses \eqn{\log \tau_{cp}} (upper) and
#' \eqn{\log \tau_{ncp}} (lower).
#'
#' Names \code{t_cp} and \code{t_ncp} kept for symmetry with the
#' Path B implementation; conceptually the returned values are
#' z-statistics under the asymptotic regime of large \code{n_boot}.
#'
#' Returns \code{NA_real_} for both statistics when the test
#' cannot be computed (any zero or negative reference scale, any
#' zero coordinate variance in the full data, length mismatch,
#' \code{n_draws} not divisible by \code{n_chains}, or non-finite
#' aggregates).
#'
#' @param effective_coef Numeric matrix of shape
#'   \code{n_draws x n_coords}.
#' @param reference_scale_per_draw Strictly positive numeric vector
#'   of length \code{n_draws}.
#' @param n_chains Integer. Number of chains in the underlying
#'   sampler.
#' @param tau_cp,tau_ncp Numeric scalars (positive).
#' @param alpha Numeric scalar. Retained for documentation; the
#'   function returns the raw z-statistics so the caller compares
#'   to \code{qnorm(1 - alpha)} explicitly.
#' @param n_boot Integer. Number of bootstrap replicates.
#' @param block_size Integer. Block length for the chain-aware
#'   bootstrap.
#' @return List with elements \code{t_cp} and \code{t_ncp}.
#' @keywords internal
#' @noRd
preflight_info_ratio_t <- function(effective_coef,
                                   reference_scale_per_draw,
                                   n_chains,
                                   tau_cp = 5, tau_ncp = 2,
                                   alpha = 0.025,
                                   n_boot = 200L,
                                   block_size = 10L) {
  if (!is.matrix(effective_coef)) {
    effective_coef <- as.matrix(effective_coef)
  }
  n_draws <- nrow(effective_coef)
  n_coords <- ncol(effective_coef)
  if (n_draws < 1L || n_coords < 1L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  if (length(reference_scale_per_draw) != n_draws) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  if (any(is.na(effective_coef))) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  if (any(is.na(reference_scale_per_draw)) ||
      any(reference_scale_per_draw <= 0)) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  full_sd <- apply(effective_coef, 2, sd)
  if (any(is.na(full_sd)) || any(full_sd <= 0)) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  full_mean_ref <- mean(reference_scale_per_draw)
  if (!is.finite(full_mean_ref) || full_mean_ref <= 0) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  log_info_full <- log(full_mean_ref / full_sd)
  m <- mean(log_info_full)
  if (!is.finite(m)) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  n_chains <- as.integer(n_chains)
  if (n_chains < 1L || n_draws %% n_chains != 0L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  block_size <- min(as.integer(block_size), n_draws %/% n_chains)
  if (block_size < 1L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  n_boot <- as.integer(n_boot)
  if (n_boot < 2L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  boot_means <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- block_bootstrap_indices(n_draws, n_chains, block_size)
    eff_b <- effective_coef[idx, , drop = FALSE]
    ref_b <- reference_scale_per_draw[idx]
    sd_b <- apply(eff_b, 2, sd)
    if (any(is.na(sd_b)) || any(sd_b <= 0)) {
      boot_means[b] <- NA_real_
      next
    }
    mean_ref_b <- mean(ref_b)
    if (!is.finite(mean_ref_b) || mean_ref_b <= 0) {
      boot_means[b] <- NA_real_
      next
    }
    boot_means[b] <- mean(log(mean_ref_b / sd_b))
  }
  boot_means <- boot_means[is.finite(boot_means)]
  if (length(boot_means) < 2L) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  se <- sd(boot_means)
  if (!is.finite(se) || se <= 0) {
    return(list(t_cp = NA_real_, t_ncp = NA_real_))
  }
  list(
    t_cp  = (m - log(tau_cp)) / se,
    t_ncp = (m - log(tau_ncp)) / se
  )
}

#' Internal: resolve parametrization flags for the long fit
#'
#' Combines the user-supplied \code{parametrization} (global) and
#' \code{parametrization_a} / \code{parametrization_W} (per-component)
#' with the optional pre-flight diagnostic to decide the CP/NCP flag
#' for each component. The pre-flight runs only if at least one
#' component still requires diagnosis after applying the user's
#' explicit overrides.
#'
#' @param parametrization Character: \code{"auto"}, \code{"ncp"} or
#'   \code{"cp"} (already validated by the caller).
#' @param parametrization_a \code{NULL} or \code{"ncp"}/\code{"cp"}.
#' @param parametrization_W \code{NULL} or \code{"ncp"}/\code{"cp"}.
#' @param prior \code{gdpar_prior} object.
#' @param stan_data List from \code{assemble_stan_data}.
#' @param amm \code{amm_spec} object.
#' @param preflight_seed Integer or \code{NULL}.
#' @param verbose Logical.
#' @return List with \code{cp_a} (logical), \code{cp_W} (logical),
#'   \code{meta} (the pre-flight metadata when run, or a stub
#'   indicating the user-explicit decision otherwise).
#' @keywords internal
#' @noRd
resolve_parametrization <- function(parametrization,
                                    parametrization_a,
                                    parametrization_W,
                                    prior, stan_data, amm,
                                    preflight_seed,
                                    verbose) {
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
    pf <- preflight_parametrization(prior, stan_data, amm,
                                    preflight_seed = preflight_seed,
                                    verbose = verbose)
    cp_a <- if (!is.na(cp_a_user)) cp_a_user else pf$cp_a
    cp_W <- if (!is.na(cp_W_user)) cp_W_user else pf$cp_W
    meta <- pf$meta
    if (!is.na(cp_a_user)) meta$decision_reason_a <- "user_explicit_a"
    if (!is.na(cp_W_user)) meta$decision_reason_W <- "user_explicit_W"
  } else {
    cp_a <- if (!is.na(cp_a_user)) cp_a_user else cp_a_global
    cp_W <- if (!is.na(cp_W_user)) cp_W_user else cp_W_global
    reason_a <- if (!is.na(cp_a_user)) "user_explicit_a"
                else "user_global"
    reason_W <- if (!is.na(cp_W_user)) "user_explicit_W"
                else "user_global"
    meta <- list(
      used_preflight       = FALSE,
      n_divergent          = NA_integer_,
      div_pct              = NA_real_,
      ebfmi_min            = NA_real_,
      t_attribution_a      = NA_real_,
      t_attribution_W      = NA_real_,
      t_info_cp_a          = NA_real_,
      t_info_ncp_a         = NA_real_,
      t_info_cp_W          = NA_real_,
      t_info_ncp_W         = NA_real_,
      decision_reason_a    = reason_a,
      decision_reason_W    = reason_W
    )
  }

  list(cp_a = as.logical(cp_a),
       cp_W = as.logical(cp_W),
       meta = meta)
}

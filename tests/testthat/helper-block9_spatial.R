# Block 9, Sub-block 9.4 / Axis 2 spatial (Session B9.37, decision D100).
#
# Frozen byte-for-byte copy of the PRE-REFACTOR monolithic
# `gdpar_dependence_robust` (Session B9.8, D75), renamed. The refactor of
# B9.37 extracted the shared `.gdpar_dependence_robust_engine`; the
# non-negotiable gate is that the temporal path stays bit-identical. The
# gated test `REFACTOR GATE: temporal path is bit-identical ...` runs this
# frozen reference against the live (engine-backed) function with the same
# seed and asserts an identical robust-SE table. It is deliberately a frozen
# copy (not a call into the live internals) so that any future drift in the
# temporal path is caught, not absorbed.
.gdpar_dependence_robust_pre_refactor <- function(object, data, index = NULL,
                                                  block_length = NULL,
                                                  type = c("moving", "circular"),
                                                  B = 199L, level = 0.95,
                                                  seed = NULL,
                                                  iter_warmup = 500L,
                                                  iter_sampling = 500L,
                                                  chains = 2L, verbose = TRUE,
                                                  ...) {
  .gdpar_assert_scalar_eb(object, "object")
  assert_data_frame(data, "data")
  type <- match.arg(type)
  assert_count(B, "B")
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_count(chains, "chains")
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort("Argument 'verbose' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  require_suggested("cmdstanr", "refit Path 1 Empirical-Bayes models")
  require_suggested("posterior", "extract posterior draws")

  n <- nrow(data)
  if (!is.null(index)) {
    if (length(index) != n) {
      gdpar_abort(
        sprintf("Argument 'index' must have length n = %d.", n),
        class = "gdpar_input_error"
      )
    }
    data <- data[order(index), , drop = FALSE]
  }
  if (is.null(block_length)) block_length <- .gdpar_default_block_length(n)
  block_length <- as.integer(block_length)
  if (block_length < 1L || block_length > n) {
    gdpar_abort(
      sprintf("Argument 'block_length' must be in [1, %d].", n),
      class = "gdpar_input_error"
    )
  }

  estimate <- .gdpar_eb_estimate_vector(object)
  model_se <- .gdpar_eb_model_se_vector(object)
  param_names <- names(estimate)
  model_se <- model_se[param_names]

  if (verbose) {
    gdpar_inform(
      sprintf(
        paste(
          "gdpar_dependence_robust is an opt-in, computationally expensive",
          "audit: it will perform %d block-bootstrap refit(s) (block length",
          "%d, %s scheme). It does not change the point estimates; it",
          "re-estimates their uncertainty to be robust to temporal",
          "dependence."
        ),
        B, block_length, type
      ),
      class = "gdpar_optin_message"
    )
  }

  if (!is.null(seed)) set.seed(as.integer(seed))
  refit_seeds <- sample.int(.Machine$integer.max, B)

  call <- object$call
  caller_env <- parent.frame()
  boot <- matrix(NA_real_, nrow = B, ncol = length(param_names),
                 dimnames = list(NULL, param_names))
  warnings_msg <- character(0L)

  for (b in seq_len(B)) {
    idx <- .gdpar_block_bootstrap_data_indices(n, block_length, type)
    sub <- data[idx, , drop = FALSE]
    sub_call <- call
    sub_call$data          <- quote(sub)
    sub_call$iter_warmup   <- iter_warmup
    sub_call$iter_sampling <- iter_sampling
    sub_call$chains        <- chains
    sub_call$verbose       <- FALSE
    sub_call$refresh       <- 0L
    sub_call$skip_id_check <- TRUE
    sub_call$seed          <- refit_seeds[b]
    env <- new.env(parent = caller_env)
    env$sub <- sub
    refit_failure_msg <- NULL
    vec <- tryCatch({
      fit_b <- eval(sub_call, envir = env)
      .gdpar_eb_estimate_vector(fit_b)
    }, error = function(e) {
      refit_failure_msg <<- sprintf(
        "Refit %d failed: %s.", b, conditionMessage(e)
      )
      NULL
    })
    if (is.null(vec)) {
      if (!is.null(refit_failure_msg)) {
        warnings_msg <- c(warnings_msg, refit_failure_msg)
      }
      next
    }
    common <- intersect(param_names, names(vec))
    boot[b, common] <- vec[common]
  }

  ok <- stats::complete.cases(boot)
  B_ok <- sum(ok)
  if (B_ok < 2L) {
    gdpar_abort(
      "Fewer than two successful bootstrap refits; cannot estimate robust SE.",
      class = "gdpar_diagnostic_error",
      data = list(B = B, B_ok = B_ok)
    )
  }
  boot_ok <- boot[ok, , drop = FALSE]
  robust_se <- apply(boot_ok, 2L, stats::sd)
  alpha <- 1 - level
  ci <- apply(boot_ok, 2L, stats::quantile,
              probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

  table <- data.frame(
    parameter = param_names,
    estimate  = unname(estimate),
    model_se  = unname(model_se),
    robust_se = unname(robust_se),
    se_ratio  = unname(robust_se / model_se),
    ci_lower  = ci[1L, ],
    ci_upper  = ci[2L, ],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  obj <- list(
    table          = table,
    block_length   = block_length,
    type           = type,
    B              = as.integer(B),
    B_ok           = as.integer(B_ok),
    level          = level,
    index_supplied = !is.null(index),
    seed           = if (is.null(seed)) NA_integer_ else as.integer(seed),
    warnings       = warnings_msg
  )
  class(obj) <- c("gdpar_dependence_robust", "list")
  obj
}

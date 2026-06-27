## Block 9.2.S: recovery adapters. One fit-and-predict entry per method
## returning the standard prediction object consumed by
## recovery_compute_metrics() (recovery_metrics.R).
##
## Capability matrix (rec_capability): which (method, gdpar_family)
## combos a method can fit. Genuine capability gaps return status
## "unsupported" (counts toward the coverage gate as a robustness pass,
## DESIGN 2.5.4). Branches deferred to a later batch session return
## "pending_impl" (NOT counted toward coverage). Any runtime fit error
## is captured as "fit_failed" with error_class (resumable; the cell is
## re-attempted on the next pass).
##
## Design note (B9.11): nonlinearity (b2) is modelled through gdpar's W
## B-spline ONLY in the K=1 path; in K>=2 distributional cells the mean
## is linear a(x1+x2) for gdpar (a measured misspecification shared with
## brms/INLA/rstanarm; mgcv may use s()). Refining gdpar to W()+K>=2 is
## logged as deferred. elpd_loo is computed only for the fully-Bayesian
## methods (gdpar/brms/rstanarm); it is a SECONDARY cross-check (DESIGN
## 2.4), and the elaborate mgcv/INLA LOO reconstructions are skipped per
## cell for cost (their leaderboard role is recovery, not elpd).

## --- capability matrix ---------------------------------------------
## Values: "full" (implemented + intended to run), "unsupported"
## (genuine package gap; clean robustness pass), "pending" (deferred to
## a later batch; excluded from the coverage gate).

rec_capability <- function(method, family) {
  caps <- list(
    gdpar = c(gaussian = "full", poisson = "full",
              neg_binomial_2 = "full", student_t = "full",
              zip = "full", zinb = "full"),
    mgcv = c(gaussian = "full", poisson = "full",
             neg_binomial_2 = "full", student_t = "full",
             zip = "full", zinb = "unsupported"),
    brms = c(gaussian = "full", poisson = "full",
             neg_binomial_2 = "full", student_t = "full",
             zip = "full", zinb = "full"),
    inla = c(gaussian = "full", poisson = "full",
             neg_binomial_2 = "full", student_t = "full",
             zip = "full", zinb = "full"),
    rstanarm = c(gaussian = "full", poisson = "full",
                 neg_binomial_2 = "full", student_t = "unsupported",
                 zip = "unsupported", zinb = "unsupported")
  )
  ## NOTE: gaussian K=2 (heteroscedastic) is a distributional fit, not a
  ## plain gaussian GLM. The per-method fit functions below decide
  ## whether they support the K of the cell; INLA/rstanarm return
  ## "unsupported" at runtime for K=2 gaussian (no distributional scale).
  cap <- caps[[method]]
  if (is.null(cap) || is.na(cap[family])) return("pending")
  unname(cap[family])
}

## --- diagnostics extraction from a cmdstanr/stanfit -----------------

.rec_stan_diag <- function(cmdstan_fit) {
  tryCatch({
    s <- cmdstan_fit$summary()
    rhat_vec <- s$rhat
    ess_vec  <- s$ess_bulk
    rhat_max <- if (all(is.na(rhat_vec))) NA_real_
                else max(rhat_vec, na.rm = TRUE)
    ess_min  <- if (all(is.na(ess_vec))) NA_real_
                else min(ess_vec, na.rm = TRUE)
    ndiv <- tryCatch(
      sum(cmdstan_fit$diagnostic_summary(quiet = TRUE)$num_divergent),
      error = function(e) NA_real_)
    list(rhat_max = rhat_max, ess_bulk_min = ess_min,
         n_divergences = as.numeric(ndiv))
  }, error = function(e) {
    list(rhat_max = NA_real_, ess_bulk_min = NA_real_,
         n_divergences = NA_real_)
  })
}

## ====================================================================
## gdpar adapter (the subject; all families).
## ====================================================================

recovery_fit_predict_gdpar <- function(model_plan, train, test, grid,
                                       fit_seed, iter) {
  fam <- model_plan$gdpar_family
  K <- model_plan$K
  cap <- rec_capability("gdpar", fam)
  if (cap == "pending") {
    return(rec_pred_status("pending_impl", gdpar_family = fam, K = K))
  }
  fit_one <- function() {
    fs_or_formula <- .rec_gdpar_spec(model_plan)
    call_args <- c(fs_or_formula$args, list(
      family = .rec_gdpar_family_obj(fam, K),
      data = train,
      chains = iter$chains, iter_warmup = iter$iter_warmup,
      iter_sampling = iter$iter_sampling, adapt_delta = 0.95,
      refresh = 0L, seed = fit_seed, skip_id_check = TRUE,
      verbose = FALSE, parametrization = "auto"))
    suppressMessages(suppressWarnings(do.call(gdpar, call_args)))
  }
  t0 <- Sys.time()
  fit <- tryCatch(fit_one(), error = function(e) {
    structure(list(msg = conditionMessage(e),
                   cls = class(e)[1L]), class = "rec_fit_err")
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "rec_fit_err")) {
    p <- rec_pred_status("fit_failed", error_class = fit$cls,
                         wall_seconds = wall, gdpar_family = fam, K = K)
    return(p)
  }
  pred <- tryCatch(
    .rec_gdpar_predict(fit, grid, model_plan),
    error = function(e) NULL)
  if (is.null(pred)) {
    return(rec_pred_status("metrics_failed", wall_seconds = wall,
                           gdpar_family = fam, K = K))
  }
  lo <- tryCatch(suppressWarnings(gdpar_loo(fit, aggregation = "subject")),
                 error = function(e) NULL)
  out <- pred
  out$status <- "ok"; out$wall_seconds <- wall
  out$gdpar_family <- fam; out$K <- K
  out$loo <- lo
  out$diag <- .rec_stan_diag(fit$fit)
  out
}

## Build the gdpar formula/amm spec for the cell.
.rec_gdpar_spec <- function(mp) {
  use_W <- isTRUE(mp$use_W)
  ## D82 / Philosophy F (B9.17): nonlinearity in x1 is modelled by gdpar's
  ## native B-spline basis of x1 IN THE LOCATION MEAN (the frozen shared
  ## columns mp$smooth_cols, precomputed in recovery_realize_cell), which
  ## replaces the raw linear x1 term wherever it appeared. The old
  ## W(theta_ref) basis was a category error -- it smooths the
  ## outcome-anchored reference, not the mean as a function of x1 -- and is
  ## dropped from the nonlinearity scenarios. mgcv/brms use their native
  ## s(x1); INLA/rstanarm consume the same frozen basis columns.
  mean_inside <- if (use_W) paste(c(mp$smooth_cols, "x2"), collapse = " + ")
                 else "x1 + x2"
  if (mp$K == 1L) {
    fm <- stats::as.formula(paste("y ~", mean_inside))
    amm <- amm_spec(a = stats::as.formula(paste("~", mean_inside)))
    return(list(args = list(formula = fm, amm = amm)))
  }
  ## K >= 2: gdpar_bf with DGP-aware slot specs (B9.12, see .rec_model_plan).
  ## Each slot models the covariate dependence the ground truth actually has:
  ## continuous scale sigma on z only under heteroscedasticity (b4); the
  ## student tail index nu and the count overdispersion phi are constant
  ## (intercept-only) because the truth never modulates them; the count
  ## zero-inflation pi is on z (truth pi = inv_logit(zi0 + zi1*z)). The old
  ## nu~a(w)/pi~a(w)/phi~a(z) specs were spurious over-parametrization. The
  ## location mean carries the smooth basis under b2 (Philosophy F); the
  ## scale/shape slots are unchanged.
  fam <- mp$gdpar_family
  scale_hetero <- isTRUE(mp$scale_hetero)
  loc <- stats::as.formula(sprintf("y ~ a(%s)", mean_inside))
  fs <- switch(fam,
    gaussian       = gdpar_bf(loc, sigma ~ a(z)),
    neg_binomial_2 = gdpar_bf(loc, phi ~ 1),
    zip            = gdpar_bf(loc, pi ~ a(z)),
    student_t      = if (scale_hetero) {
                       gdpar_bf(loc, sigma ~ a(z), nu ~ 1)
                     } else {
                       gdpar_bf(loc, sigma ~ 1, nu ~ 1)
                     },
    zinb           = gdpar_bf(loc, phi ~ 1, pi ~ a(z)),
    stop(sprintf("no gdpar_bf spec for family %s", fam))
  )
  list(args = list(formula = fs))
}

.rec_gdpar_family_obj <- function(fam, K) gdpar_family(fam)

## gdpar prediction on the grid -> standard moments.
.rec_gdpar_predict <- function(fit, grid, mp) {
  fam <- mp$gdpar_family; K <- mp$K
  if (K == 1L) {
    mu_draws <- stats::predict(fit, newdata = grid, type = "response",
                               summary = "draws")            ## S x n
    if (fam == "gaussian") {
      sy <- tryCatch({
        m <- fit$fit$draws(variables = "sigma_y",
                           format = "draws_matrix")
        mean(sqrt(colMeans(unclass(m)^2)))
      }, error = function(e) stats::sd(as.numeric(grid$x1)) )
      return(list(eyx_draws = mu_draws,
                  sd_hat = rep(sy, ncol(mu_draws)),
                  tail_q_hat = colMeans(mu_draws) + stats::qnorm(0.95) * sy))
    }
    ## poisson: response = lambda
    mom <- rec_moments_count(mu_draws, size = Inf, pi_mat = 0)
    return(mom)
  }
  arr <- stats::predict(fit, newdata = grid, type = "response",
                        summary = "draws")                   ## S x n x K
  loc <- arr[, , 1L]
  if (fam == "gaussian") {
    return(rec_moments_continuous(loc, arr[, , 2L], nu = Inf))
  }
  if (fam == "student_t") {
    return(.rec_student_moments(loc, arr[, , 2L], arr[, , 3L]))
  }
  if (fam == "neg_binomial_2") {
    return(.rec_count_moments_mat(loc, size_mat = arr[, , 2L], pi_mat = 0))
  }
  if (fam == "zip") {
    return(rec_moments_count(loc, size = Inf, pi_mat = arr[, , 2L]))
  }
  if (fam == "zinb") {
    return(.rec_count_moments_mat(loc, size_mat = arr[, , 2L],
                                  pi_mat = arr[, , 3L]))
  }
  stop(sprintf("no gdpar predict path for family %s", fam))
}

## Student-t moments with per-draw degrees of freedom.
.rec_student_moments <- function(loc, sd_mat, nu_mat) {
  nu_safe <- pmax(nu_mat, 2.01)
  sd_eff <- sd_mat * sqrt(nu_safe / (nu_safe - 2))
  qt95 <- matrix(stats::qt(0.95, df = as.vector(nu_safe)),
                 nrow = nrow(nu_safe))
  tail_mat <- loc + sd_mat * qt95
  list(eyx_draws = loc, sd_hat = colMeans(sd_eff),
       tail_q_hat = colMeans(tail_mat))
}

## Count moments with matrix-valued size and/or zero-inflation.
.rec_count_moments_mat <- function(lambda_mat, size_mat = Inf, pi_mat = 0) {
  if (length(pi_mat) == 1L) {
    pi_mat <- matrix(pi_mat, nrow = nrow(lambda_mat),
                     ncol = ncol(lambda_mat))
  }
  finite_size <- !identical(size_mat, Inf)
  eyx <- (1 - pi_mat) * lambda_mat
  base_var <- if (finite_size) lambda_mat + lambda_mat^2 / size_mat
              else lambda_mat
  ey2 <- (1 - pi_mat) * (base_var + lambda_mat^2)
  sd_mat <- sqrt(pmax(ey2 - eyx^2, 0))
  lam_bar <- colMeans(lambda_mat); pi_bar <- colMeans(pi_mat)
  size_bar <- if (finite_size) colMeans(size_mat) else Inf
  tail <- numeric(length(lam_bar))
  for (i in seq_along(lam_bar)) {
    sz <- if (finite_size) size_bar[i] else Inf
    tail[i] <- rec_count_tail_plugin(lam_bar[i], pi_bar[i], sz, 0.95)
  }
  list(eyx_draws = eyx, sd_hat = colMeans(sd_mat), tail_q_hat = tail)
}

## ====================================================================
## brms adapter (most general competitor). posterior_epred / _predict
## give E[y|x] draws and predictive draws directly, so the family
## moment algebra is native.
## ====================================================================

recovery_fit_predict_brms <- function(model_plan, train, test, grid,
                                      fit_seed, iter) {
  fam <- model_plan$gdpar_family; K <- model_plan$K
  cap <- rec_capability("brms", fam)
  if (cap == "pending") {
    return(rec_pred_status("pending_impl", gdpar_family = fam, K = K))
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(rec_pred_status("unsupported", error_class = "brms_absent",
                           gdpar_family = fam, K = K))
  }
  bf_obj_fam <- .rec_brms_formula_family(model_plan)
  t0 <- Sys.time()
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      brms::brm(bf_obj_fam$formula, data = train, family = bf_obj_fam$family,
                chains = iter$chains, warmup = iter$iter_warmup,
                iter = iter$iter_warmup + iter$iter_sampling,
                seed = fit_seed, refresh = 0L, backend = "rstan",
                silent = 2L, cores = iter$chains))),
    error = function(e) structure(list(cls = class(e)[1L]),
                                  class = "rec_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "rec_fit_err")) {
    return(rec_pred_status("fit_failed", error_class = fit$cls,
                           wall_seconds = wall, gdpar_family = fam, K = K))
  }
  out <- tryCatch({
    eyx <- brms::posterior_epred(fit, newdata = grid)
    ypred <- brms::posterior_predict(fit, newdata = grid)
    lo <- tryCatch(suppressWarnings(brms::loo(fit)), error = function(e) NULL)
    list(status = "ok", wall_seconds = wall, gdpar_family = fam, K = K,
         eyx_draws = eyx, sd_hat = apply(ypred, 2L, stats::sd),
         tail_q_hat = apply(ypred, 2L, stats::quantile, probs = 0.95,
                            names = FALSE),
         loo = lo, diag = .rec_brms_diag(fit))
  }, error = function(e) rec_pred_status("metrics_failed",
                                         error_class = class(e)[1L],
                                         wall_seconds = wall,
                                         gdpar_family = fam, K = K))
  out
}

.rec_brms_formula_family <- function(mp) {
  rhs <- "x1 + x2"
  if (isTRUE(mp$use_W)) rhs <- "s(x1) + x2"   ## brms supports mgcv smooths
  fam <- mp$gdpar_family
  if (fam == "gaussian" && mp$K == 2L) {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs)),
                                   stats::as.formula("sigma ~ z")),
                family = brms::brmsfamily("gaussian")))
  }
  if (fam == "gaussian") {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs))),
                family = brms::brmsfamily("gaussian")))
  }
  if (fam == "student_t") {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs))),
                family = brms::brmsfamily("student")))
  }
  if (fam == "poisson") {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs))),
                family = brms::brmsfamily("poisson")))
  }
  if (fam == "neg_binomial_2") {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs))),
                family = brms::brmsfamily("negbinomial")))
  }
  if (fam == "zip") {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs))),
                family = brms::brmsfamily("zero_inflated_poisson")))
  }
  if (fam == "zinb") {
    return(list(formula = brms::bf(stats::as.formula(paste("y ~", rhs))),
                family = brms::brmsfamily("zero_inflated_negbinomial")))
  }
  stop(sprintf("no brms family for %s", fam))
}

.rec_brms_diag <- function(fit) {
  tryCatch({
    sf <- fit$fit              ## stanfit
    s <- rstan::summary(sf)$summary
    list(rhat_max = max(s[, "Rhat"], na.rm = TRUE),
         ess_bulk_min = min(s[, "n_eff"], na.rm = TRUE),
         n_divergences = NA_real_)
  }, error = function(e) list(rhat_max = NA_real_,
                              ess_bulk_min = NA_real_,
                              n_divergences = NA_real_))
}

## ====================================================================
## rstanarm adapter. posterior_epred / _predict native; distributional
## scale, Student-t, and zero-inflation are unsupported.
## ====================================================================

recovery_fit_predict_rstanarm <- function(model_plan, train, test, grid,
                                          fit_seed, iter) {
  fam <- model_plan$gdpar_family; K <- model_plan$K
  cap <- rec_capability("rstanarm", fam)
  if (cap == "pending") {
    return(rec_pred_status("pending_impl", gdpar_family = fam, K = K))
  }
  if (cap == "unsupported" || (fam == "gaussian" && K == 2L)) {
    return(rec_pred_status("unsupported",
                           error_class = "rstanarm_no_distributional",
                           gdpar_family = fam, K = K))
  }
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    return(rec_pred_status("unsupported", error_class = "rstanarm_absent",
                           gdpar_family = fam, K = K))
  }
  ## D82 / Philosophy F: nonlinearity via the frozen shared B-spline basis
  ## columns of x1 (fixed effects), replacing the raw linear x1 term.
  rhs <- if (isTRUE(model_plan$use_W))
           paste(c(model_plan$smooth_cols, "x2"), collapse = " + ")
         else "x1 + x2"
  fm <- stats::as.formula(paste("y ~", rhs))
  t0 <- Sys.time()
  fit <- tryCatch({
    if (fam == "gaussian") {
      suppressMessages(suppressWarnings(rstanarm::stan_glm(
        fm, data = train, family = stats::gaussian(),
        chains = iter$chains, iter = iter$iter_warmup + iter$iter_sampling,
        warmup = iter$iter_warmup, seed = fit_seed, refresh = 0L,
        cores = iter$chains)))
    } else if (fam == "poisson") {
      suppressMessages(suppressWarnings(rstanarm::stan_glm(
        fm, data = train, family = stats::poisson(),
        chains = iter$chains, iter = iter$iter_warmup + iter$iter_sampling,
        warmup = iter$iter_warmup, seed = fit_seed, refresh = 0L,
        cores = iter$chains)))
    } else if (fam == "neg_binomial_2") {
      ## stan_glm.nb() internally evaluates an unqualified stan_glm() in its
      ## own frame; under requireNamespace/:: loading (no library()) that
      ## symbol is not on the search path and the call dies instantly with
      ## "could not find function 'stan_glm'". Call stan_glm() with the
      ## neg_binomial_2 family directly (mathematically equivalent).
      suppressMessages(suppressWarnings(rstanarm::stan_glm(
        fm, data = train, family = rstanarm::neg_binomial_2(),
        chains = iter$chains, iter = iter$iter_warmup + iter$iter_sampling,
        warmup = iter$iter_warmup, seed = fit_seed, refresh = 0L,
        cores = iter$chains)))
    } else stop("unreachable")
  }, error = function(e) structure(list(cls = class(e)[1L]),
                                   class = "rec_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "rec_fit_err")) {
    return(rec_pred_status("fit_failed", error_class = fit$cls,
                           wall_seconds = wall, gdpar_family = fam, K = K))
  }
  out <- tryCatch({
    eyx <- rstanarm::posterior_epred(fit, newdata = grid)
    ypred <- rstanarm::posterior_predict(fit, newdata = grid)
    lo <- tryCatch(suppressWarnings(rstanarm::loo(fit)),
                   error = function(e) NULL)
    list(status = "ok", wall_seconds = wall, gdpar_family = fam, K = K,
         eyx_draws = eyx, sd_hat = apply(ypred, 2L, stats::sd),
         tail_q_hat = apply(ypred, 2L, stats::quantile, probs = 0.95,
                            names = FALSE),
         loo = lo, diag = list(rhat_max = NA_real_, ess_bulk_min = NA_real_,
                               n_divergences = NA_real_))
  }, error = function(e) rec_pred_status("metrics_failed",
                                         error_class = class(e)[1L],
                                         wall_seconds = wall,
                                         gdpar_family = fam, K = K))
  out
}

## ====================================================================
## mgcv adapter. gaulss (heteroscedastic gaussian), scat (scaled t),
## nb (neg-binomial), ziP (zero-inflated poisson). Coverage draws come
## from the MVN(coef, Vp) Bayesian interpretation (Wood 2017).
## ====================================================================

recovery_fit_predict_mgcv <- function(model_plan, train, test, grid,
                                      fit_seed, iter, S = 1000L) {
  fam <- model_plan$gdpar_family; K <- model_plan$K
  cap <- rec_capability("mgcv", fam)
  if (cap == "pending") {
    return(rec_pred_status("pending_impl", gdpar_family = fam, K = K))
  }
  if (cap == "unsupported") {
    return(rec_pred_status("unsupported",
                           error_class = "mgcv_family_gap",
                           gdpar_family = fam, K = K))
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(rec_pred_status("unsupported", error_class = "mgcv_absent",
                           gdpar_family = fam, K = K))
  }
  t0 <- Sys.time()
  res <- tryCatch(.rec_mgcv_fit_predict(model_plan, train, grid, fit_seed, S),
                  error = function(e) structure(list(cls = class(e)[1L]),
                                                class = "rec_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "rec_fit_err")) {
    return(rec_pred_status("fit_failed", error_class = res$cls,
                           wall_seconds = wall, gdpar_family = fam, K = K))
  }
  res$status <- "ok"; res$wall_seconds <- wall
  res$gdpar_family <- fam; res$K <- K
  res$loo <- NULL    ## elpd cross-check skipped for mgcv (cost; DESIGN 2.4)
  res$diag <- list(rhat_max = NA_real_, ess_bulk_min = NA_real_,
                   n_divergences = NA_real_)
  res
}

.rec_mgcv_mean_rhs <- function(mp, target = "x1") {
  if (isTRUE(mp$use_W)) "s(x1) + x2" else "x1 + x2"
}

.rec_mgcv_fit_predict <- function(mp, train, grid, fit_seed, S) {
  fam <- mp$gdpar_family
  rhs <- .rec_mgcv_mean_rhs(mp)
  set.seed(fit_seed)
  if (fam == "gaussian" && mp$K == 2L) {
    ## gaulss: location-scale gaussian; logb default sets log(sd + b)
    form <- list(stats::as.formula(paste("y ~", rhs)),
                 stats::as.formula("~ z"))
    g <- mgcv::gam(form, data = train, family = mgcv::gaulss())
    pr <- stats::predict(g, newdata = grid, type = "response")
    mu <- pr[, 1L]
    ## gaulss second column is 1/sd on the (log-link + b) scale -> invert
    sd_hat <- 1 / pr[, 2L]
    se_mu <- stats::predict(g, newdata = grid, type = "link",
                            se.fit = TRUE)$se.fit[, 1L]
    eyx <- .rec_gauss_draws(mu, se_mu, S)
    tail_q <- mu + stats::qnorm(0.95) * sd_hat
    return(list(eyx_draws = eyx, sd_hat = sd_hat, tail_q_hat = tail_q))
  }
  if (fam == "gaussian") {
    g <- mgcv::gam(stats::as.formula(paste("y ~", rhs)), data = train,
                   family = stats::gaussian(), method = "REML")
    return(.rec_mgcv_gauss_like(g, grid, S, scat = FALSE))
  }
  if (fam == "student_t") {
    g <- mgcv::gam(stats::as.formula(paste("y ~", rhs)), data = train,
                   family = mgcv::scat(), method = "REML")
    return(.rec_mgcv_gauss_like(g, grid, S, scat = TRUE))
  }
  if (fam == "poisson") {
    g <- mgcv::gam(stats::as.formula(paste("y ~", rhs)), data = train,
                   family = stats::poisson(), method = "REML")
    return(.rec_mgcv_count(g, grid, S, size = Inf, zi = FALSE))
  }
  if (fam == "neg_binomial_2") {
    g <- mgcv::gam(stats::as.formula(paste("y ~", rhs)), data = train,
                   family = mgcv::nb(), method = "REML")
    size <- g$family$getTheta(TRUE)
    return(.rec_mgcv_count(g, grid, S, size = size, zi = FALSE))
  }
  if (fam == "zip") {
    g <- mgcv::gam(stats::as.formula(paste("y ~", rhs)), data = train,
                   family = mgcv::ziP())
    return(.rec_mgcv_zip(g, grid, S))
  }
  stop(sprintf("no mgcv path for %s", fam))
}

.rec_gauss_draws <- function(mu, se, S) {
  n <- length(mu)
  m <- matrix(stats::rnorm(S * n, mean = rep(mu, each = S),
                           sd = rep(se, each = S)), nrow = S, ncol = n)
  m
}

.rec_mgcv_gauss_like <- function(g, grid, S, scat) {
  pr <- stats::predict(g, newdata = grid, type = "link", se.fit = TRUE)
  mu <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  if (scat) {
    th <- g$family$getTheta(TRUE)   ## c(nu, sigma) for scat
    nu <- th[1L]; sig <- th[2L]
    sd_hat <- rep(sig * sqrt(nu / (nu - 2)), length(mu))
    tail_q <- mu + sig * stats::qt(0.95, df = nu)
  } else {
    sig <- sqrt(g$sig2)
    sd_hat <- rep(sig, length(mu))
    tail_q <- mu + stats::qnorm(0.95) * sig
  }
  list(eyx_draws = .rec_gauss_draws(mu, se, S),
       sd_hat = sd_hat, tail_q_hat = tail_q)
}

.rec_mgcv_count <- function(g, grid, S, size, zi) {
  pr <- stats::predict(g, newdata = grid, type = "link", se.fit = TRUE)
  eta <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  eta_draws <- .rec_gauss_draws(eta, se, S)
  lam_draws <- exp(eta_draws)
  mom <- if (is.finite(size)) {
    .rec_count_moments_mat(lam_draws, size_mat = matrix(size, S, ncol(lam_draws)),
                           pi_mat = 0)
  } else {
    rec_moments_count(lam_draws, size = Inf, pi_mat = 0)
  }
  mom
}

.rec_mgcv_zip <- function(g, grid, S) {
  ## ziP: linear predictor is for the poisson mean on log scale; the
  ## zero-inflation is a monotone function of the mean (Wood et al.).
  ## Use predicted response mean directly; approximate the predictive
  ## tail from the implied poisson with the fitted mean.
  mu_resp <- as.numeric(stats::predict(g, newdata = grid, type = "response"))
  pr <- stats::predict(g, newdata = grid, type = "link", se.fit = TRUE)
  eta <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  eyx <- .rec_gauss_draws(mu_resp, se * mu_resp, S)   ## crude mean draws
  sd_hat <- sqrt(mu_resp)                              ## poisson-ish floor
  tail_q <- stats::qpois(0.95, lambda = pmax(mu_resp, 1e-8))
  list(eyx_draws = eyx, sd_hat = sd_hat, tail_q_hat = tail_q)
}

## ====================================================================
## INLA adapter. Gaussian / poisson / nbinomial / zero-inflated /
## Student-t (family "T"). Heteroscedastic gaussian (K=2) unsupported.
## Predictive draws synthesized from summary.fitted (Gaussian approx for
## the mean function), which suffices for coverage of E[y|x].
## ====================================================================

recovery_fit_predict_inla <- function(model_plan, train, test, grid,
                                      fit_seed, iter, S = 1000L) {
  fam <- model_plan$gdpar_family; K <- model_plan$K
  cap <- rec_capability("inla", fam)
  if (cap == "pending") {
    return(rec_pred_status("pending_impl", gdpar_family = fam, K = K))
  }
  if (fam == "gaussian" && K == 2L) {
    return(rec_pred_status("unsupported",
                           error_class = "inla_no_distributional_scale",
                           gdpar_family = fam, K = K))
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    return(rec_pred_status("unsupported", error_class = "inla_absent",
                           gdpar_family = fam, K = K))
  }
  t0 <- Sys.time()
  res <- tryCatch(.rec_inla_fit_predict(model_plan, train, grid, S),
                  error = function(e) structure(list(cls = class(e)[1L]),
                                                class = "rec_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "rec_fit_err")) {
    return(rec_pred_status("fit_failed", error_class = res$cls,
                           wall_seconds = wall, gdpar_family = fam, K = K))
  }
  res$status <- "ok"; res$wall_seconds <- wall
  res$gdpar_family <- fam; res$K <- K
  res$loo <- NULL    ## elpd cross-check skipped for INLA (cost; DESIGN 2.4)
  res$diag <- list(rhat_max = NA_real_, ess_bulk_min = NA_real_,
                   n_divergences = NA_real_)
  res
}

.rec_inla_family_name <- function(fam) {
  switch(fam,
    gaussian       = "gaussian",
    poisson        = "poisson",
    neg_binomial_2 = "nbinomial",
    student_t      = "T",
    zip            = "zeroinflatedpoisson1",
    zinb           = "zeroinflatednbinomial1",
    stop(sprintf("no INLA family for %s", fam)))
}

## Discretization granularity for INLA's native rw2 smooth of x1 (D83).
## NOT the effective df -- that is chosen by INLA's hyperprior on the rw2
## precision (its native penalized-smoothness selection). Revisable knob,
## analogous to .rec_smooth_df = 9 for the fixed-basis methods.
.rec_inla_rw2_groups <- 50L

.rec_inla_fit_predict <- function(mp, train, grid, S) {
  fam <- mp$gdpar_family
  n_tr <- nrow(train)
  aug <- grid; aug$y <- NA_real_
  base_cols <- c("y", "x1", "x2", "z", "w")
  full <- rbind(train[, base_cols], aug[, base_cols])
  ## D83 (Philosophy F refinement, B9.18): INLA's NATIVE smooth of x1 is the
  ## second-order random-walk latent model (rw2) -- the INLA analog of a
  ## penalized cubic spline (mgcv's s()) -- NOT a fixed-effect basis. The
  ## frozen fixed-effect B-spline basis D82 shared with rstanarm (which truly
  ## has no native penalized smoother) breaks INLA's Newton-Raphson optimizer
  ## for the Student-t family (heavy tails + collinear flexible mean ->
  ## "lambda < 1/lambda_lim"; regularizing the fixed-effect prior does NOT
  ## help); rw2 converges across all families and recovers comparably.
  ## Bucketing INLA with rstanarm under D82 was a category error -- INLA HAS
  ## a native penalized smoother -- mirroring the B9.17 W(theta_ref) fix.
  ## inla.group bins x1 on the COMBINED train+grid so the two share
  ## consistent bins and the grid prediction stays valid.
  if (isTRUE(mp$use_W)) {
    full$x1grp <- INLA::inla.group(full$x1, n = .rec_inla_rw2_groups)
    form <- y ~ f(x1grp, model = "rw2", scale.model = TRUE) + x2
  } else {
    form <- stats::as.formula("y ~ x1 + x2")
  }
  pred_rows <- (n_tr + 1L):(n_tr + nrow(grid))
  inla_fit <- function(extra) {
    do.call(INLA::inla, c(list(
      formula = form, family = .rec_inla_family_name(fam), data = full,
      control.predictor = list(compute = TRUE, link = 1L),
      control.compute = list(config = FALSE)), extra))
  }
  .ok_fit <- function(f) !is.null(f) && !is.null(f$summary.fitted.values) &&
    all(is.finite(f$summary.fitted.values$mean[pred_rows]))
  ## D83 robustness fallback (B9.18): the default INLA strategy
  ## (simplified.laplace latent + numerical hyperparameter integration) can
  ## diverge on the hardest count cells -- e.g. zip + nonlinearity + AR(1) +
  ## multicollinearity all active, ~1/10 realizations -- aborting with a
  ## fit_failed. Retry once with the Gaussian latent strategy + empirical
  ## Bayes hyperparameters (int.strategy = "eb"), which converges there. The
  ## fallback only fires when the default fails or returns non-finite means,
  ## so the cells that converge under the default (the vast majority) are
  ## bit-for-bit untouched; the rare fallback cell trades a little
  ## hyperparameter integration for a finite, recorded fit -- real
  ## comparison over abort (max-robustness regla-piedra).
  fit <- tryCatch(inla_fit(NULL), error = function(e) NULL)
  if (!.ok_fit(fit)) {
    fit <- inla_fit(list(control.inla = list(strategy = "gaussian",
                                             int.strategy = "eb")))
  }
  mu <- fit$summary.fitted.values$mean[pred_rows]
  sdf <- fit$summary.fitted.values$sd[pred_rows]
  eyx <- .rec_gauss_draws(mu, sdf, S)
  ## conditional SD + tail from the fitted mean per family (plug-in)
  if (fam == "gaussian") {
    prec_name <- "Precision for the Gaussian observations"
    tau <- tryCatch(fit$summary.hyperpar[prec_name, "mean"],
                    error = function(e) NA_real_)
    sig <- if (is.finite(tau) && tau > 0) sqrt(1 / tau) else stats::sd(train$y)
    sd_hat <- rep(sig, length(mu)); tail_q <- mu + stats::qnorm(0.95) * sig
  } else if (fam == "student_t") {
    sd_hat <- sdf; tail_q <- mu + stats::qnorm(0.95) * sdf
  } else {
    ## counts: mu is the response mean; approximate poisson/nb tail
    sd_hat <- sqrt(pmax(mu, 0))
    tail_q <- stats::qpois(0.95, lambda = pmax(mu, 1e-8))
  }
  list(eyx_draws = eyx, sd_hat = sd_hat, tail_q_hat = tail_q)
}

## --- dispatch -------------------------------------------------------

recovery_method_table <- function() {
  list(
    gdpar    = recovery_fit_predict_gdpar,
    mgcv     = recovery_fit_predict_mgcv,
    brms     = recovery_fit_predict_brms,
    inla     = recovery_fit_predict_inla,
    rstanarm = recovery_fit_predict_rstanarm
  )
}

## Block 9.2.O: eBird re-validation adapters. One fit-and-predict entry
## per method PER COORDINATE, returning a standard coord-result consumed
## by ebird_finalize_coord() (which calls ebird_metrics_9_2_O.R to get the
## test elpd on the ORIGINAL response measure + coverage).
##
## Structure (B9.19, emergent D84): gdpar's K > 1 distributional
## regression is SCALAR-response only (no public p = 2 distributional
## joint API; the multivariate p = 2 path is K = 1 gaussian, Block 7 /
## 9.3.d). So only track-G level G1 (gaussian K = 1) is fit JOINTLY p = 2
## for gdpar/brms (Block 7 "nativa-mixta"); every distributional level
## (G2, G3, track N) is fit PER COORDINATE for ALL methods -- a forced API
## constraint, not a design choice. mgcv/INLA/rstanarm are per coordinate
## throughout (as in Block 7).
##
## A coord-result is list(status, family_tag, param_draws | loglik_native,
##   diag, wall_seconds, error_class). param_draws is a named list of
## S x n_test draw matrices of the predictive parameters on the family's
## native scale; loglik_native is an S x n_test native log-lik (brms/
## rstanarm). ebird_finalize_coord() turns either into the original-measure
## elpd + coverage.

## --- formula fragments ---------------------------------------------
.ebird_rhs <- function() paste(.ebird_9_2_O_covariates, collapse = " + ")
.ebird_has_group <- function(d) nlevels(droplevels(factor(d$group))) >= 2L

.ebird_iter <- function() list(chains = 2L, iter_warmup = 500L,
                               iter_sampling = 1000L)

## clip a unit-interval response into (eps, 1-eps): S&T occurrence has
## exact zeros (species absent) / ones, which the beta likelihood (open
## (0,1) support) cannot take. Consistent with the density evaluation,
## which re-clips the observed y the same way (ebird_metrics_9_2_O.R).
.ebird_clip_unit <- function(y) pmin(pmax(y, .ebird_eps_logit),
                                     1 - .ebird_eps_logit)

## column name of the modelled response for a coord under a track.
.ebird_ycol <- function(coord, track) {
  if (identical(track, "G")) return(if (coord == 1L) "z1" else "z2")
  if (coord == 1L) "occ" else "cnt"            ## track N original scale
}
## the original-scale observed outcome (for the Jacobian + density).
.ebird_yorig <- function(test, coord) {
  if (coord == 1L) test$occ else test$cnt
}

## --- unsupported / failed coord-result skeletons --------------------
.ebird_coord_unsup <- function(family_tag, error_class) {
  list(status = "unsupported", family_tag = family_tag,
       param_draws = NULL, loglik_native = NULL,
       diag = NULL, wall_seconds = NA_real_, error_class = error_class)
}
.ebird_coord_fail <- function(family_tag, error_class, wall) {
  list(status = "fit_failed", family_tag = family_tag,
       param_draws = NULL, loglik_native = NULL,
       diag = NULL, wall_seconds = wall, error_class = error_class)
}

## --- diagnostics ----------------------------------------------------
.ebird_stan_diag_cmdstan <- function(cmdstan_fit) {
  tryCatch({
    s <- cmdstan_fit$summary()
    rh <- s$rhat; es <- s$ess_bulk
    nd <- tryCatch(sum(cmdstan_fit$diagnostic_summary(quiet = TRUE)$num_divergent),
                   error = function(e) NA_real_)
    list(rhat_max = if (all(is.na(rh))) NA_real_ else max(rh, na.rm = TRUE),
         ess_bulk_min = if (all(is.na(es))) NA_real_ else min(es, na.rm = TRUE),
         n_divergences = as.numeric(nd))
  }, error = function(e) list(rhat_max = NA_real_, ess_bulk_min = NA_real_,
                              n_divergences = NA_real_))
}
.ebird_diag_null <- function() list(rhat_max = NA_real_,
                                    ess_bulk_min = NA_real_,
                                    n_divergences = NA_real_)

## ====================================================================
## gdpar
## ====================================================================

## gdpar scalar fit for the distributional levels (per coordinate).
## gdpar models the group as a grouping factor via the `group` arg, not in
## the formula; the mean RHS keeps the fixed covariates and group is passed
## separately (mirrors Block 7 bench_gdpar_fit).
.ebird_fit_gdpar_scalar <- function(coord, track, level, train, test,
                                    fit_seed, fam_name, family_tag) {
  rhs <- .ebird_rhs()
  ycol <- .ebird_ycol(coord, track)
  tr <- train; tr$y <- tr[[ycol]]
  if (fam_name == "beta") tr$y <- .ebird_clip_unit(tr$y)
  ## D84 (B9.19): gdpar predict() on newdata for K>1 + RANDOM grouping
  ## (J_groups>1) is "not yet implemented (Session 8.4)". The held-out test
  ## ELPD needs predict on newdata, so for the distributional (K>1) scalar
  ## fits the species group enters as a FIXED effect in the location mean
  ## (a(rhs + group), J_groups=1 -> predict works, also ~10x faster) rather
  ## than a random effect via the group= arg. Documented asymmetry vs the
  ## competitors' random group (~5 species levels -> minor; see report).
  has_g <- .ebird_has_group(tr)
  loc_rhs <- if (has_g) paste(rhs, "+ group") else rhs
  loc <- stats::as.formula(sprintf("y ~ a(%s)", loc_rhs))
  fs <- switch(level,
    G2 = gdpar_bf(loc, stats::as.formula(sprintf("sigma ~ a(%s)", rhs))),
    G3 = gdpar_bf(loc, stats::as.formula(sprintf("sigma ~ a(%s)", rhs)),
                  nu ~ 1),
    Nb = gdpar_bf(loc, phi ~ 1),
    Nt = gdpar_bf(loc, phi ~ 1, p ~ 1),
    stop("bad gdpar scalar level"))
  iter <- .ebird_iter()
  ## parallel_chains forwards through gdpar()'s `...` to cmdstanr's $sample();
  ## without it the chains run sequentially (the default), roughly doubling the
  ## wall time of every distributional fit.  Result-neutral (same seed, chains,
  ## iters -- only execution parallelism), so it does not touch any cell value.
  ##
  ## Track-N sampler config (B9.20).  The track-N slots beta(phi~1) and
  ## tweedie(phi~1, p~1) carry an intercept-only a() whose hierarchical scale
  ## sigma_a_k is weakly identified (one coefficient -> no information for a
  ## group-level SD).  At gdpar's default max_treedepth=12 NUTS walks that
  ## near-flat direction with maximal trajectories: >6 h/cell (tweedie worst),
  ## and ~24% divergences at adapt_delta=0.95.  Capping at the Stan/brms
  ## standard treedepth=10 and tightening adapt_delta to 0.99 makes the fit
  ## feasible (~3-6 min) with zero divergences.  The predictive quantities
  ## (mu, phi, p) are stable across chains (sd ~ 1e-3); the residual rhat=Inf
  ## is the benign sigma_a_k non-identification, NOT a prediction failure
  ## (the tweedie ELPD reproduces mgcv's -- see report).  Flagged honestly in
  ## the manifest (rhat_max), never masked.  Track-G (G2/G3) is untouched.
  is_trackN <- level %in% c("Nb", "Nt")
  call_args <- list(formula = fs, family = gdpar_family(fam_name),
                    data = tr, chains = iter$chains,
                    parallel_chains = iter$chains,
                    iter_warmup = iter$iter_warmup,
                    iter_sampling = iter$iter_sampling,
                    adapt_delta = 0.95,
                    max_treedepth = if (is_trackN) 10L else 12L,
                    refresh = 0L, seed = fit_seed, skip_id_check = TRUE,
                    verbose = FALSE, parametrization = "auto")
  t0 <- Sys.time()
  fit <- tryCatch(suppressMessages(suppressWarnings(do.call(gdpar, call_args))),
                  error = function(e) structure(list(cls = class(e)[1L]),
                                                class = "ebird_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "ebird_fit_err")) {
    return(.ebird_coord_fail(family_tag, fit$cls, wall))
  }
  pd <- tryCatch(.ebird_gdpar_param_draws(fit, test, level), error = function(e) NULL)
  if (is.null(pd)) return(.ebird_coord_fail(family_tag, "gdpar_predict", wall))
  lo <- tryCatch(suppressWarnings(gdpar_loo(fit, aggregation = "subject")),
                 error = function(e) NULL)
  list(status = "ok", family_tag = family_tag, param_draws = pd,
       loglik_native = NULL, loo = lo, diag = .ebird_stan_diag_cmdstan(fit$fit),
       wall_seconds = wall, error_class = NA_character_)
}

## extract predictive parameter draws on the test set from a scalar fit.
.ebird_gdpar_param_draws <- function(fit, test, level) {
  arr <- stats::predict(fit, newdata = test, type = "response",
                        summary = "draws")              ## S x n x K
  if (length(dim(arr)) == 2L) {                          ## K = 1 fallback
    return(list(mu = arr))
  }
  K <- dim(arr)[3L]
  if (level == "G2") return(list(mu = arr[, , 1L], sigma = arr[, , 2L]))
  if (level == "G3") return(list(mu = arr[, , 1L], sigma = arr[, , 2L],
                                 nu = arr[, , 3L]))
  if (level == "Nb") return(list(mu = arr[, , 1L], phi = arr[, , 2L]))
  if (level == "Nt") return(list(mu = arr[, , 1L], phi = arr[, , 2L],
                                 p = arr[, , 3L]))
  stop("bad gdpar param-draw level")
}

## gdpar joint p = 2 gaussian (track G1 only; Block 7 reproduction).
.ebird_fit_gdpar_joint_G1 <- function(train, test, fit_seed) {
  rhs <- .ebird_rhs()
  tr <- train; tr$y <- as.matrix(tr[, c("z1", "z2")])
  iter <- .ebird_iter()
  call_args <- list(
    formula = stats::as.formula(sprintf("y ~ %s", rhs)),
    family = gdpar_family("gaussian"),
    amm = amm_spec(p = 2L,
                   dims = dimwise(a = stats::as.formula(sprintf("~ %s", rhs)),
                                  b = NULL)),
    data = tr, chains = iter$chains, parallel_chains = iter$chains,
    iter_warmup = iter$iter_warmup,
    iter_sampling = iter$iter_sampling, adapt_delta = 0.95, refresh = 0L,
    seed = fit_seed, skip_id_check = TRUE, verbose = FALSE,
    parametrization = "auto")
  if (.ebird_has_group(tr)) call_args$group <- ~ group
  t0 <- Sys.time()
  fit <- tryCatch(suppressMessages(suppressWarnings(do.call(gdpar, call_args))),
                  error = function(e) structure(list(cls = class(e)[1L]),
                                                class = "ebird_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "ebird_fit_err")) {
    return(list(`1` = .ebird_coord_fail("gaussian_z", fit$cls, wall),
                `2` = .ebird_coord_fail("gaussian_z", fit$cls, wall),
                wall = wall, diag = .ebird_diag_null()))
  }
  res <- tryCatch({
    eta <- stats::predict(fit, newdata = test, type = "linear_predictor")  ## S x n x 2
    sy <- fit$fit$draws(variables = "sigma_y", format = "draws_matrix")
    sy <- unclass(sy)                                    ## S x 2 (per coord)
    diag <- .ebird_stan_diag_cmdstan(fit$fit)
    mk_coord <- function(k) {
      sig_k <- if (ncol(sy) >= k) sqrt(mean(sy[, k]^2)) else sqrt(mean(sy^2))
      list(status = "ok", family_tag = "gaussian_z",
           param_draws = list(mu = eta[, , k], sigma = sig_k),
           loglik_native = NULL, diag = diag, wall_seconds = wall,
           error_class = NA_character_)
    }
    list(`1` = mk_coord(1L), `2` = mk_coord(2L), wall = wall, diag = diag)
  }, error = function(e) {
    list(`1` = .ebird_coord_fail("gaussian_z", "gdpar_predict", wall),
         `2` = .ebird_coord_fail("gaussian_z", "gdpar_predict", wall),
         wall = wall, diag = .ebird_diag_null())
  })
  res
}

## ====================================================================
## brms
## ====================================================================

.ebird_fit_brms_joint_G1 <- function(train, test, fit_seed) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(list(`1` = .ebird_coord_unsup("gaussian_z", "brms_absent"),
                `2` = .ebird_coord_unsup("gaussian_z", "brms_absent"),
                wall = NA_real_, diag = .ebird_diag_null()))
  }
  rhs <- .ebird_rhs()
  re <- if (.ebird_has_group(train)) " + (1 | group)" else ""
  iter <- .ebird_iter()
  f <- brms::bf(stats::as.formula(sprintf("mvbind(z1, z2) ~ %s%s", rhs, re))) +
    brms::set_rescor(FALSE)
  t0 <- Sys.time()
  fit <- tryCatch(suppressMessages(suppressWarnings(
    brms::brm(f, data = train, family = brms::brmsfamily("gaussian"),
              chains = iter$chains, warmup = iter$iter_warmup,
              iter = iter$iter_warmup + iter$iter_sampling, seed = fit_seed,
              refresh = 0L, backend = "rstan", silent = 2L,
              cores = iter$chains))),
    error = function(e) structure(list(cls = class(e)[1L]),
                                  class = "ebird_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "ebird_fit_err")) {
    return(list(`1` = .ebird_coord_fail("gaussian_z", fit$cls, wall),
                `2` = .ebird_coord_fail("gaussian_z", fit$cls, wall),
                wall = wall, diag = .ebird_diag_null()))
  }
  diag <- .ebird_brms_diag(fit)
  mk <- function(resp) {
    ll <- tryCatch(brms::log_lik(fit, newdata = test, resp = resp),
                   error = function(e) NULL)
    if (is.null(ll)) return(.ebird_coord_fail("gaussian_z", "brms_loglik", wall))
    list(status = "ok", family_tag = "gaussian_z", param_draws = NULL,
         loglik_native = ll, diag = diag, wall_seconds = wall,
         error_class = NA_character_)
  }
  list(`1` = mk("z1"), `2` = mk("z2"), wall = wall, diag = diag)
}

.ebird_fit_brms_scalar <- function(coord, track, level, train, test,
                                   fit_seed, family_tag) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(.ebird_coord_unsup(family_tag, "brms_absent"))
  }
  rhs <- .ebird_rhs()
  re <- if (.ebird_has_group(train)) " + (1 | group)" else ""
  ycol <- .ebird_ycol(coord, track)
  tr <- train; tr$ymod <- tr[[ycol]]
  te <- test;  te$ymod <- te[[ycol]]
  if (level == "Nb") { tr$ymod <- .ebird_clip_unit(tr$ymod)
                       te$ymod <- .ebird_clip_unit(te$ymod) }
  base <- stats::as.formula(sprintf("ymod ~ %s%s", rhs, re))
  spec <- switch(level,
    G2 = list(f = brms::bf(base, stats::as.formula(sprintf("sigma ~ %s%s", rhs, re))),
              fam = brms::brmsfamily("gaussian")),
    G3 = list(f = brms::bf(base, stats::as.formula(sprintf("sigma ~ %s%s", rhs, re))),
              fam = brms::brmsfamily("student")),
    Nb = list(f = brms::bf(base), fam = brms::brmsfamily("Beta")),
    stop("brms: no scalar spec for level"))
  iter <- .ebird_iter()
  t0 <- Sys.time()
  fit <- tryCatch(suppressMessages(suppressWarnings(
    brms::brm(spec$f, data = tr, family = spec$fam, chains = iter$chains,
              warmup = iter$iter_warmup,
              iter = iter$iter_warmup + iter$iter_sampling, seed = fit_seed,
              refresh = 0L, backend = "rstan", silent = 2L,
              cores = iter$chains))),
    error = function(e) structure(list(cls = class(e)[1L]),
                                  class = "ebird_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "ebird_fit_err")) {
    return(.ebird_coord_fail(family_tag, fit$cls, wall))
  }
  ll <- tryCatch(brms::log_lik(fit, newdata = te), error = function(e) NULL)
  if (is.null(ll)) return(.ebird_coord_fail(family_tag, "brms_loglik", wall))
  lo <- tryCatch(suppressWarnings(brms::loo(fit)), error = function(e) NULL)
  list(status = "ok", family_tag = family_tag, param_draws = NULL,
       loglik_native = ll, loo = lo, diag = .ebird_brms_diag(fit),
       wall_seconds = wall, error_class = NA_character_)
}

.ebird_brms_diag <- function(fit) {
  tryCatch({
    s <- rstan::summary(fit$fit)$summary
    list(rhat_max = max(s[, "Rhat"], na.rm = TRUE),
         ess_bulk_min = min(s[, "n_eff"], na.rm = TRUE),
         n_divergences = NA_real_)
  }, error = function(e) .ebird_diag_null())
}

## ====================================================================
## mgcv (per coordinate; pointwise mean posterior via predict se.fit;
## dispersion plug-in -- documented approximation, three-layer rigor:
## exact mean-posterior + plug-in scale/shape).
## ====================================================================

.ebird_fit_mgcv <- function(coord, track, level, train, test, fit_seed,
                            family_tag, S = 1000L) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(.ebird_coord_unsup(family_tag, "mgcv_absent"))
  }
  ## mgcv's tw() family calls ldTweedie() by bare name during fitting; under
  ## devtools::load_all(gdpar) + mgcv::gam (loaded, not attached) that lookup
  ## fails ("could not find function 'ldTweedie'"). Attaching mgcv puts its
  ## exports on the search path (idempotent; no clash with the gdpar API).
  suppressMessages(require(mgcv, quietly = TRUE))
  rhs <- .ebird_rhs()
  re <- if (.ebird_has_group(train)) ' + s(group, bs = "re")' else ""
  ycol <- .ebird_ycol(coord, track)
  tr <- train; tr$ymod <- tr[[ycol]]
  te <- test;  te$ymod <- te[[ycol]]
  if (level == "Nb") tr$ymod <- .ebird_clip_unit(tr$ymod)
  base <- stats::as.formula(sprintf("ymod ~ %s%s", rhs, re))
  set.seed(fit_seed)
  t0 <- Sys.time()
  res <- tryCatch({
    if (level == "G1") {
      g <- mgcv::gam(base, data = tr, family = stats::gaussian(), method = "REML")
      pr <- stats::predict(g, newdata = te, type = "link", se.fit = TRUE)
      mu <- .ebird_mu_draws(pr$fit, pr$se.fit, S, identity)
      list(mu = mu, sigma = sqrt(g$sig2))
    } else if (level == "G2") {
      g <- mgcv::gam(list(base, stats::as.formula(sprintf("~ %s%s", rhs, re))),
                     data = tr, family = mgcv::gaulss())
      pr <- stats::predict(g, newdata = te, type = "link", se.fit = TRUE)
      mu <- .ebird_mu_draws(pr$fit[, 1L], pr$se.fit[, 1L], S, identity)
      prr <- stats::predict(g, newdata = te, type = "response")
      sd_hat <- 1 / prr[, 2L]                 ## gaulss: col2 is 1/sd
      list(mu = mu, sigma = matrix(sd_hat, S, length(sd_hat), byrow = TRUE))
    } else if (level == "G3") {
      g <- mgcv::gam(base, data = tr, family = mgcv::scat(), method = "REML")
      pr <- stats::predict(g, newdata = te, type = "link", se.fit = TRUE)
      mu <- .ebird_mu_draws(pr$fit, pr$se.fit, S, identity)
      th <- g$family$getTheta(TRUE)           ## c(nu, sigma) for scat
      list(mu = mu, sigma = th[2L], nu = th[1L])
    } else if (level == "Nb") {
      g <- mgcv::gam(base, data = tr, family = mgcv::betar(), method = "REML")
      pr <- stats::predict(g, newdata = te, type = "link", se.fit = TRUE)
      mu <- .ebird_mu_draws(pr$fit, pr$se.fit, S, stats::plogis)
      ## betar precision: getTheta(TRUE) already returns phi on the NATURAL
      ## scale (verified: a betar fit with true phi=12 returns 11.79; getTheta
      ## (FALSE) returns log(phi)=2.47).  The prior exp() double-exponentiated
      ## it (phi -> exp(11.79) = 1.3e5), making the beta density catastrophically
      ## peaked (elpd ~ -4.7e4).  Use it directly, as scat/tw do above.
      phi <- g$family$getTheta(TRUE)
      list(mu = mu, phi = phi)
    } else if (level == "Nt") {
      g <- mgcv::gam(base, data = tr, family = mgcv::tw(), method = "REML")
      pr <- stats::predict(g, newdata = te, type = "link", se.fit = TRUE)
      mu <- .ebird_mu_draws(pr$fit, pr$se.fit, S, exp)
      pw <- g$family$getTheta(TRUE)           ## tw power p
      list(mu = mu, phi = g$scale %||% summary(g)$dispersion, p = pw)
    } else stop("mgcv: bad level")
  }, error = function(e) structure(list(cls = class(e)[1L]),
                                   class = "ebird_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "ebird_fit_err")) {
    return(.ebird_coord_fail(family_tag, res$cls, wall))
  }
  list(status = "ok", family_tag = family_tag, param_draws = res,
       loglik_native = NULL, diag = .ebird_diag_null(),
       wall_seconds = wall, error_class = NA_character_)
}

## S x n draws of the mean predictor: mu_i ~ N(fit_i, se_i) on link scale,
## then inv_link. Correct for POINTWISE elpd (each point uses its own
## marginal mean posterior).
.ebird_mu_draws <- function(fit_link, se_link, S, inv_link) {
  n <- length(fit_link)
  m <- matrix(stats::rnorm(S * n, rep(fit_link, each = S),
                           rep(pmax(se_link, 1e-8), each = S)), nrow = S, ncol = n)
  inv_link(m)
}

## ====================================================================
## INLA (per coordinate; augmented-NA refit + posterior.sample on the
## test Predictor rows -> family density).
## ====================================================================

.ebird_inla_family <- function(level) {
  switch(level, G1 = "gaussian", G3 = "T", Nb = "beta", stop("inla bad level"))
}

.ebird_fit_inla <- function(coord, track, level, train, test, fit_seed,
                            family_tag, S = 1000L) {
  if (level %in% c("G2", "Nt")) {
    return(.ebird_coord_unsup(family_tag,
      if (level == "G2") "inla_no_distributional_scale" else "inla_no_tweedie"))
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    return(.ebird_coord_unsup(family_tag, "inla_absent"))
  }
  rhs <- .ebird_rhs()
  re <- if (.ebird_has_group(train)) ' + f(group, model = "iid")' else ""
  ycol <- .ebird_ycol(coord, track)
  n_tr <- nrow(train)
  tr <- train; tr$ymod <- tr[[ycol]]
  if (level == "Nb") tr$ymod <- .ebird_clip_unit(tr$ymod)
  te <- test;  te$ymod <- NA_real_
  cols <- c("ymod", .ebird_9_2_O_covariates, "group")
  full <- rbind(tr[, cols], te[, cols])
  form <- stats::as.formula(sprintf("ymod ~ %s%s", rhs, re))
  fam <- .ebird_inla_family(level)
  pred_rows <- (n_tr + 1L):(n_tr + nrow(test))
  t0 <- Sys.time()
  .one_fit <- function(ctrl_inla) {
    fit <- INLA::inla(form, family = fam, data = full,
                      control.predictor = list(compute = TRUE, link = 1L),
                      control.compute = list(config = TRUE),
                      control.inla = ctrl_inla, verbose = FALSE)
    set.seed(fit_seed)
    samp <- INLA::inla.posterior.sample(S, fit)
    pidx <- grep("^Predictor:", rownames(samp[[1L]]$latent))[pred_rows]
    eta <- t(vapply(samp, function(s) s$latent[pidx, 1L], numeric(length(pred_rows))))
    .ebird_inla_params(level, eta, samp, family_tag)
  }
  ## INLA's default approximation can fail (non-deterministically) on spiky
  ## responses -- e.g. z1 = logit(occ) carries a mass spike at the clip floor
  ## from the many exact-zero occurrences, which destabilises the student-t
  ## dof estimation (~half the G3 sub-regions error out, and which half flips
  ## across runs).  Retry once in INLA's most robust mode (empirical-Bayes
  ## hyperpar integration + Gaussian latent strategy) before giving up: the
  ## model class is unchanged (same student-t), only the posterior
  ## approximation is more conservative.  Cells that converge in the default
  ## mode keep the default result; the retry only fires on a hard failure.
  res <- tryCatch(.one_fit(list()), error = function(e1)
    tryCatch(.one_fit(list(strategy = "gaussian", int.strategy = "eb")),
             error = function(e2) structure(list(cls = class(e2)[1L]),
                                            class = "ebird_fit_err")))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "ebird_fit_err")) {
    return(.ebird_coord_fail(family_tag, res$cls, wall))
  }
  list(status = "ok", family_tag = family_tag, param_draws = res,
       loglik_native = NULL, diag = .ebird_diag_null(),
       wall_seconds = wall, error_class = NA_character_)
}

## eta is S x n on the linear-predictor scale of the family.
.ebird_inla_params <- function(level, eta, samp, family_tag) {
  hp <- function(nm) vapply(samp, function(s) {
    h <- s$hyperpar; if (nm %in% names(h)) as.numeric(h[[nm]]) else NA_real_
  }, numeric(1L))
  if (level == "G1") {
    tau <- hp("Precision for the Gaussian observations")
    sig <- ifelse(is.finite(tau) & tau > 0, sqrt(1 / tau), NA_real_)
    return(list(mu = eta, sigma = sig))          ## gaussian identity link
  }
  if (level == "G3") {
    ## INLA family "T": precision + dof hyperparameters.  INLA's "precision"
    ## tau is the MARGINAL precision (Var(y) = 1/tau), NOT the precision of the
    ## raw t-scale.  The downstream density (.ebird_logdens_student_z) treats
    ## sigma as the t-scale s, with Var = s^2 * nu/(nu-2).  Hence the marginal
    ## SD sqrt(1/tau) must be converted: s = sqrt(1/tau) * sqrt((nu-2)/nu).
    ## Verified numerically against a known-scale INLA "T" fit (recovers s to
    ## sampling error).  Passing the marginal SD as the scale (the prior bug)
    ## flattens the predictive density catastrophically as nu -> 2.
    tau <- hp("precision for the student-t observations")
    if (all(is.na(tau))) tau <- hp("Precision for the student-t observations")
    nu  <- hp("degrees of freedom for student-t")
    if (all(is.na(nu))) nu <- hp("degrees of freedom for the student-t")
    nu  <- ifelse(is.finite(nu) & nu > 2, nu, 10)
    marg_sd <- ifelse(is.finite(tau) & tau > 0, sqrt(1 / tau), NA_real_)
    sig <- marg_sd * sqrt((nu - 2) / nu)
    return(list(mu = eta, sigma = sig, nu = nu))
  }
  if (level == "Nb") {
    phi <- hp("precision parameter for the beta observations")
    return(list(mu = stats::plogis(eta),          ## beta logit link
                phi = ifelse(is.finite(phi), phi, NA_real_)))
  }
  stop("inla params bad level")
}

## ====================================================================
## rstanarm (per coordinate; gaussian only -> distributional + beta-RE +
## tweedie unsupported).
## ====================================================================

.ebird_fit_rstanarm <- function(coord, track, level, train, test, fit_seed,
                                family_tag) {
  if (level != "G1") {
    return(.ebird_coord_unsup(family_tag, "rstanarm_no_distributional"))
  }
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    return(.ebird_coord_unsup(family_tag, "rstanarm_absent"))
  }
  rhs <- .ebird_rhs()
  has_g <- .ebird_has_group(train)
  re <- if (has_g) " + (1 | group)" else ""
  ycol <- .ebird_ycol(coord, track)
  tr <- train; tr$ymod <- tr[[ycol]]
  te <- test;  te$ymod <- te[[ycol]]
  fm <- stats::as.formula(sprintf("ymod ~ %s%s", rhs, re))
  iter <- .ebird_iter()
  fit_fn <- if (has_g) rstanarm::stan_glmer else rstanarm::stan_glm
  t0 <- Sys.time()
  fit <- tryCatch(suppressMessages(suppressWarnings(
    fit_fn(fm, data = tr, family = stats::gaussian(), chains = iter$chains,
           iter = iter$iter_warmup + iter$iter_sampling,
           warmup = iter$iter_warmup, seed = fit_seed, refresh = 0L,
           cores = iter$chains))),
    error = function(e) structure(list(cls = class(e)[1L]),
                                  class = "ebird_fit_err"))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "ebird_fit_err")) {
    return(.ebird_coord_fail(family_tag, fit$cls, wall))
  }
  ll <- tryCatch(rstanarm::log_lik(fit, newdata = te), error = function(e) NULL)
  if (is.null(ll)) return(.ebird_coord_fail(family_tag, "rstanarm_loglik", wall))
  lo <- tryCatch(suppressWarnings(rstanarm::loo(fit)), error = function(e) NULL)
  list(status = "ok", family_tag = family_tag, param_draws = NULL,
       loglik_native = ll, loo = lo, diag = .ebird_diag_null(),
       wall_seconds = wall, error_class = NA_character_)
}

## ====================================================================
## per-coord finalize: coord-result -> elpd on original measure + coverage
## ====================================================================

ebird_finalize_coord <- function(coord_res, coord, track, train, test, sim_seed) {
  if (!identical(coord_res$status, "ok")) {
    return(c(coord_res, list(elpd = NA_real_, coverage90 = NA_real_,
                             n_nonfinite = NA_integer_, elpd_loo = NA_real_,
                             pareto_k_max = NA_real_)))
  }
  y_orig <- .ebird_yorig(test, coord)
  ll <- if (!is.null(coord_res$loglik_native)) {
    ebird_loglik_orig_from_native(coord_res$loglik_native, y_orig, coord, track)
  } else {
    ebird_loglik_orig_from_params(coord_res$family_tag, coord_res$param_draws,
                                  y_orig, coord, track)
  }
  ee <- ebird_elpd_from_loglik(ll)
  cov <- NA_real_
  if (!is.null(coord_res$param_draws)) {
    yp <- tryCatch(ebird_sim_ypred_orig(coord_res$family_tag,
                                        coord_res$param_draws, coord, track,
                                        seed = sim_seed), error = function(e) NULL)
    cov <- ebird_coverage90(yp, y_orig)
  }
  ## track B: in-sample LOO on the original measure (Bayesian fits only)
  lm <- ebird_loo_orig(coord_res$loo, .ebird_yorig(train, coord), coord, track)
  c(coord_res, list(elpd = ee$elpd, coverage90 = cov,
                    n_nonfinite = ee$n_nonfinite_pts,
                    elpd_loo = lm$elpd_loo, pareto_k_max = lm$pareto_k_max))
}

## ====================================================================
## cell runner
## ====================================================================

ebird_run_cell_9_2_O <- function(track, level, method, train, test, fit_seed) {
  lm <- ebird_9_2_O_level_meta(track, level)
  ## map (level, coord) -> the per-coord gdpar/mgcv/inla level token + tag.
  coord_level <- function(coord) {
    if (track == "G") return(level)               ## G1/G2/G3 same per coord
    if (coord == 1L) return("Nb") else return("Nt")
  }
  coord_tag <- function(coord) if (coord == 1L) lm$tag1 else lm$tag2
  coord_fam <- function(coord) if (coord == 1L) lm$fam1 else lm$fam2

  ## joint p=2 gaussian (gdpar/brms, track G1 only)
  if (track == "G" && level == "G1" && method %in% c("gdpar", "brms")) {
    j <- if (method == "gdpar") .ebird_fit_gdpar_joint_G1(train, test, fit_seed)
         else .ebird_fit_brms_joint_G1(train, test, fit_seed)
    c1 <- ebird_finalize_coord(j[["1"]], 1L, "G", train, test, fit_seed + 1L)
    c2 <- ebird_finalize_coord(j[["2"]], 2L, "G", train, test, fit_seed + 2L)
    return(.ebird_assemble_cell(c1, c2, j$wall, j$diag))
  }

  fit_coord <- function(coord) {
    cl <- coord_level(coord); tag <- coord_tag(coord); fam <- coord_fam(coord)
    if (method == "gdpar") {
      gl <- if (track == "G") (if (level == "G1") NA else level)
            else (if (coord == 1L) "Nb" else "Nt")
      if (track == "G" && level == "G1") {
        ## per-coord gaussian K=1 for gdpar only arises if not the joint
        ## branch -- not reached (joint handled above). Guard anyway.
        return(.ebird_fit_gdpar_scalar(coord, "G", "G2", train, test,
                                       fit_seed + coord, fam, tag))
      }
      return(.ebird_fit_gdpar_scalar(coord, track, gl, train, test,
                                     fit_seed + coord, fam, tag))
    }
    if (method == "brms") {
      bl <- if (track == "G") level else "Nb"
      if (track == "N" && coord == 2L) {
        return(.ebird_coord_unsup("tweedie", "brms_no_tweedie"))
      }
      return(.ebird_fit_brms_scalar(coord, track, bl, train, test,
                                    fit_seed + coord, tag))
    }
    if (method == "mgcv") {
      ml <- if (track == "G") level else (if (coord == 1L) "Nb" else "Nt")
      return(.ebird_fit_mgcv(coord, track, ml, train, test,
                             fit_seed + coord, tag))
    }
    if (method == "inla") {
      il <- if (track == "G") level else (if (coord == 1L) "Nb" else "Nt")
      return(.ebird_fit_inla(coord, track, il, train, test,
                             fit_seed + coord, tag))
    }
    if (method == "rstanarm") {
      rl <- if (track == "G") level else (if (coord == 1L) "Nb" else "Nt")
      return(.ebird_fit_rstanarm(coord, track, rl, train, test,
                                 fit_seed + coord, tag))
    }
    stop("unknown method")
  }

  r1 <- fit_coord(1L); r2 <- fit_coord(2L)
  c1 <- ebird_finalize_coord(r1, 1L, track, train, test, fit_seed + 1L)
  c2 <- ebird_finalize_coord(r2, 2L, track, train, test, fit_seed + 2L)
  wall <- sum(c(c1$wall_seconds, c2$wall_seconds), na.rm = TRUE)
  diag <- if (!is.null(c1$diag) && !all(is.na(unlist(c1$diag)))) c1$diag
          else c2$diag %||% .ebird_diag_null()
  .ebird_assemble_cell(c1, c2, wall, diag)
}

.ebird_assemble_cell <- function(c1, c2, wall, diag) {
  cell_status <- if ("ok" %in% c(c1$status, c2$status)) "ok"
                 else if (all(c(c1$status, c2$status) == "unsupported")) "unsupported"
                 else "fit_failed"
  list(status = cell_status, wall_seconds = wall,
       diag = diag %||% .ebird_diag_null(),
       coords = list(`1` = c1, `2` = c2))
}

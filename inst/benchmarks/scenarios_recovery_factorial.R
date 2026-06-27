## Block 9.2.S: synthetic recovery factorial with KNOWN ground truth.
##
## Emits the complete 2 x 2^5 = 64 scenario factorial of Decision D'
## (canonized B9.10, D78). Each scenario fixes a data-generating
## process (DGP) with a-priori known structural parameters, so the
## recovery metrics (recovery_metrics.R) can compare each method's
## estimated mean / scale / tail against the TRUE conditional moments
## on a fixed evaluation grid.
##
## Two response regimes, each a 2^5 factorial over five pathologies
## (one on/off bit each):
##   regime C (continuous): b1 multicollinearity, b2 nonlinearity,
##     b3 autocorrelation, b4 heteroscedasticity, b5 heavy tails.
##   regime K (count):      b1 multicollinearity, b2 nonlinearity,
##     b3 autocorrelation, b4 zero-inflation, b5 overdispersion.
## The three covariate-side bits (b1,b2,b3) are shared; the two
## response-side bits (b4,b5) are regime-specific (DESIGN_9_2_EXTERNAL
## section 2). bitmask string = "b1b2b3b4b5"; scenario_id = "C_10010".
##
## Honesty (no overreach, feedback_no_overreach_uniqueness): gdpar does
## NOT model autocorrelation. In b3-active cells the recovery measures
## ROBUSTNESS to misspecification (finite, non-collapsing recovery),
## not recovery of a modelled structure; mgcv/INLA may model it and are
## expected to dominate there. The autocorrelation enters as additive
## AR(1) error (continuous) or an AR(1) latent log-rate perturbation
## (count); the TRUE mean function is the structural (AR-free) mean.
##
## Determinism (DESIGN section 3.4; NOT the frozen golden band
## 91001-91099, per D76): truth seed per scenario in 91005..91068;
## per-cell realization and grid seeds DERIVED from it. Statistical
## quantities over R reps are reproducible but not frozen bit-exact.

## ---------------------------------------------------------------------
## Fixed structural constants (the ground truth). Identical across reps;
## only the data realization (noise + covariate draw) varies per rep.
## ---------------------------------------------------------------------

.rec_const <- list(
  ## continuous mean: mu(x) = b0 + b1*x1 + b2*x2 ( + nonlin(x1) )
  cont_b0 = 0.30, cont_b1 = 0.80, cont_b2 = -0.50,
  ## count log-mean: log lambda(x) = c0 + c1*x1 + c2*x2 ( + nonlin(x1) )
  count_c0 = 0.50, count_c1 = 0.45, count_c2 = -0.30,
  ## nonlinearity term added to the (continuous) mean / (count) log-mean
  nonlin_amp = 0.90, nonlin_freq = 1.6,
  count_nonlin_amp = 0.45, count_nonlin_freq = 1.4,
  ## homoscedastic continuous SD (b4 off)
  sigma0 = 0.50,
  ## heteroscedastic continuous log-SD: log sigma(z) = ls0 + ls1*z
  ls0 = log(0.45), ls1_on = 0.55, ls1_off = 0.0,
  ## heavy-tail Student-t degrees of freedom (b5 on)
  nu_heavy = 3.0,
  ## autocorrelation AR(1) coefficient (b3 on)
  ar_rho = 0.60,
  ## multicollinearity: x2 = x1 + N(0, mc_sd^2) when b1 on
  mc_sd = 0.05,
  ## count overdispersion neg-binomial size (b5 on): smaller = more disp
  nb_size = 2.0,
  ## count zero-inflation logit: pi(z) = inv_logit(zi0 + zi1*z) (b4 on)
  zi0 = -0.60, zi1 = 0.50,
  ## evaluation-grid size and base seed band
  n_grid = 200L,
  truth_seed_base = 91005L
)

.rec_inv_logit <- function(x) 1 / (1 + exp(-x))

## Per-scenario truth seed: contiguous 91005..91068. Index 0..63 with
## continuous regime first (0..31), then count (32..63).
.rec_truth_seed <- function(scenario_index) {
  .rec_const$truth_seed_base + as.integer(scenario_index)
}

## Derived seeds (DESIGN 3.4). grid fixed per scenario; realization per
## (scenario, rep); HMC fit seed per (scenario, rep, method).
.rec_grid_seed   <- function(truth_seed) truth_seed * 1000L + 999L
.rec_realiz_seed <- function(truth_seed, rep) truth_seed * 1000L + rep * 10L
rec_fit_seed     <- function(truth_seed, rep, method_idx) {
  truth_seed * 1000L + rep * 10L + as.integer(method_idx)
}

## Method index map (used for fit-seed derivation; stable order).
.rec_method_idx <- c(gdpar = 1L, mgcv = 2L, brms = 3L,
                     inla = 4L, rstanarm = 5L)

## ---------------------------------------------------------------------
## Scenario table: 64 rows with id / regime / bitmask / pathologies /
## truth_seed / gdpar_family / K. The gdpar family + K are determined by
## the active response-side bits (b4,b5) per regime.
## ---------------------------------------------------------------------

recovery_scenario_table <- function() {
  bits_grid <- expand.grid(b5 = c(0L, 1L), b4 = c(0L, 1L),
                           b3 = c(0L, 1L), b2 = c(0L, 1L),
                           b1 = c(0L, 1L))
  bits_grid <- bits_grid[, c("b1", "b2", "b3", "b4", "b5")]
  bits_grid <- bits_grid[order(bits_grid$b1, bits_grid$b2, bits_grid$b3,
                               bits_grid$b4, bits_grid$b5), ]
  rows <- list()
  idx <- 0L
  for (regime in c("C", "K")) {
    for (r in seq_len(nrow(bits_grid))) {
      b <- bits_grid[r, ]
      bitmask <- paste0(b$b1, b$b2, b$b3, b$b4, b$b5)
      fam_K <- .rec_family_for(regime, b$b4, b$b5)
      paths <- .rec_pathologies_active(regime, b)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario_index    = idx,
        scenario_id       = sprintf("%s_%s", regime, bitmask),
        regime            = regime,
        bitmask           = bitmask,
        b1 = b$b1, b2 = b$b2, b3 = b$b3, b4 = b$b4, b5 = b$b5,
        pathologies_active = paths,
        truth_seed        = .rec_truth_seed(idx),
        gdpar_family      = fam_K$family,
        K                 = fam_K$K,
        stringsAsFactors  = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

.rec_family_for <- function(regime, b4, b5) {
  if (regime == "C") {
    if (b5 == 1L) return(list(family = "student_t", K = 3L))
    if (b4 == 1L) return(list(family = "gaussian",  K = 2L))
    return(list(family = "gaussian", K = 1L))
  }
  ## count regime
  if (b4 == 1L && b5 == 1L) return(list(family = "zinb",            K = 3L))
  if (b4 == 1L)             return(list(family = "zip",             K = 2L))
  if (b5 == 1L)             return(list(family = "neg_binomial_2",  K = 2L))
  list(family = "poisson", K = 1L)
}

.rec_pathologies_active <- function(regime, b) {
  labs_cov <- c(multicollinearity = b$b1, nonlinearity = b$b2,
                autocorrelation = b$b3)
  labs_resp <- if (regime == "C") {
    c(heteroscedasticity = b$b4, heavy_tails = b$b5)
  } else {
    c(zero_inflation = b$b4, overdispersion = b$b5)
  }
  active <- c(labs_cov, labs_resp)
  on <- names(active)[active == 1L]
  if (length(on) == 0L) "none" else paste(on, collapse = "+")
}

## ---------------------------------------------------------------------
## Covariate draw (shared by realization and grid). Columns x1, x2, z, w.
## Multicollinearity (b1) makes x2 nearly aliased with x1.
## ---------------------------------------------------------------------

.rec_draw_covariates <- function(n, b1, seed) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- if (b1 == 1L) x1 + stats::rnorm(n, sd = .rec_const$mc_sd)
        else stats::rnorm(n)
  z  <- stats::rnorm(n)
  w  <- stats::rnorm(n)
  data.frame(x1 = x1, x2 = x2, z = z, w = w)
}

## Structural linear predictor of the LOCATION slot (mean for continuous,
## log-rate for count), AR-free. Used both to realize data and to define
## the true mean on the grid.
.rec_eta_location <- function(cov, regime, b2) {
  C <- .rec_const
  if (regime == "C") {
    eta <- C$cont_b0 + C$cont_b1 * cov$x1 + C$cont_b2 * cov$x2
    if (b2 == 1L) eta <- eta + C$nonlin_amp * sin(C$nonlin_freq * cov$x1)
  } else {
    eta <- C$count_c0 + C$count_c1 * cov$x1 + C$count_c2 * cov$x2
    if (b2 == 1L) eta <- eta +
      C$count_nonlin_amp * sin(C$count_nonlin_freq * cov$x1)
  }
  eta
}

## Structural SD of the continuous response (AR-free, marginal).
.rec_sigma_cont <- function(cov, b4) {
  C <- .rec_const
  if (b4 == 1L) exp(C$ls0 + C$ls1_on * cov$z) else rep(C$sigma0, nrow(cov))
}

## Structural zero-inflation probability (count, b4).
.rec_pi_count <- function(cov, b4) {
  C <- .rec_const
  if (b4 == 1L) .rec_inv_logit(C$zi0 + C$zi1 * cov$z)
  else rep(0, nrow(cov))
}

## ---------------------------------------------------------------------
## AR(1) error / latent generator (b3). Returns a length-n AR(1) series
## with unit marginal variance (continuous) scaled by the caller, or a
## modest-variance latent (count).
## ---------------------------------------------------------------------

.rec_ar1_series <- function(n, rho, innov_sd) {
  e <- numeric(n)
  e[1L] <- stats::rnorm(1L, sd = innov_sd / sqrt(1 - rho^2))
  for (t in 2:n) e[t] <- rho * e[t - 1L] + stats::rnorm(1L, sd = innov_sd)
  e
}

## ---------------------------------------------------------------------
## Nonlinearity (b2) shared smooth basis (B9.17, decision D82, Philosophy
## F "comparacion real"; refined by D83, B9.18). The DGP adds
## amp*sin(freq*x1) to the (log-)mean. Each method recovers it with its
## NATIVE smooth of x1: mgcv/brms keep their penalized s(x1); INLA uses its
## native second-order random walk f(x1, model="rw2") (D83; see the INLA
## adapter); gdpar uses a(bs(x1)) in the mean (its native model.matrix
## path); rstanarm -- the only method with no native penalized smoother --
## consumes a FROZEN, SHARED B-spline basis of x1 as plain fixed-effect
## columns. The basis is fitted on the training x1 once and applied
## identically to test/grid (predict on the bs object), so train and the
## evaluation grid share the exact same knots -- this sidesteps every
## fixed-basis method's predict-time basis-reconstruction quirk and gives
## gdpar and rstanarm an identical basis (they then differ only by their
## regularization). df = 9 (degree-3 B-spline) matches mgcv's default
## k = 10 basis dimension; revisable knob. D83 corrected D82's category
## error of bucketing INLA with rstanarm: INLA HAS a native penalized
## smoother (rw2), and the frozen fixed-effect basis broke its
## Newton-Raphson optimizer for the Student-t family. The even older
## W(theta_ref) basis was itself a category error (smooth on the
## outcome-anchored reference, not the mean-in-x1) and was dropped in B9.17.
## ---------------------------------------------------------------------

.rec_smooth_df     <- 9L
.rec_smooth_degree <- 3L
.rec_smooth_cols   <- function(df = .rec_smooth_df) paste0("x1_s", seq_len(df))

## Augment train/test/grid with the frozen B-spline basis columns of x1.
.rec_augment_bspline <- function(train, test, grid,
                                 df = .rec_smooth_df,
                                 degree = .rec_smooth_degree) {
  bs_obj <- suppressWarnings(
    splines::bs(train$x1, df = df, degree = degree))
  cols <- .rec_smooth_cols(df)
  apply_basis <- function(d) {
    B <- suppressWarnings(stats::predict(bs_obj, d$x1))
    B <- matrix(as.numeric(B), nrow = nrow(d))
    colnames(B) <- cols
    cbind(d, as.data.frame(B))
  }
  list(train = apply_basis(train), test = apply_basis(test),
       grid = apply_basis(grid))
}

## ---------------------------------------------------------------------
## Realize one cell's data + grid + ground-truth moments + model plan.
## scenario_row: one row of recovery_scenario_table(). rep: 1..R.
## Returns: train / test / grid (data.frames), truth_grid (mu, sd,
## tail_q095 vectors on grid), model_plan (per-method modeling recipe),
## scenario_meta.
## ---------------------------------------------------------------------

recovery_realize_cell <- function(scenario_row, rep, n = 2000L,
                                  n_train = NULL) {
  stopifnot(nrow(scenario_row) == 1L)
  ## 80/20 train/test split scaled to n (2000 -> 1600/400, DESIGN 2.3)
  if (is.null(n_train)) n_train <- as.integer(round(0.8 * n))
  regime  <- scenario_row$regime
  b1 <- scenario_row$b1; b2 <- scenario_row$b2; b3 <- scenario_row$b3
  b4 <- scenario_row$b4; b5 <- scenario_row$b5
  truth_seed <- scenario_row$truth_seed
  C <- .rec_const

  ## --- realize observed data (per rep) ---
  realiz_seed <- .rec_realiz_seed(truth_seed, rep)
  cov <- .rec_draw_covariates(n, b1, realiz_seed)
  set.seed(realiz_seed + 7L)  ## separate stream for the noise
  eta_loc <- .rec_eta_location(cov, regime, b2)
  y <- .rec_simulate_response(regime, eta_loc, cov, b3, b4, b5)
  full <- data.frame(y = y, x1 = cov$x1, x2 = cov$x2, z = cov$z, w = cov$w)
  ## also expose y1 for adapters that index outcomes by coordinate
  full$y1 <- full$y
  train <- full[seq_len(n_train), , drop = FALSE]
  test  <- full[(n_train + 1L):n, , drop = FALSE]

  ## --- evaluation grid + ground-truth moments (per scenario) ---
  grid_seed <- .rec_grid_seed(truth_seed)
  grid_cov <- .rec_draw_covariates(C$n_grid, b1, grid_seed)
  grid <- data.frame(x1 = grid_cov$x1, x2 = grid_cov$x2,
                     z = grid_cov$z, w = grid_cov$w)
  truth_grid <- .rec_truth_moments(regime, grid_cov, b2, b4, b5)

  ## Nonlinearity (b2): augment with the frozen shared B-spline basis of
  ## x1 (D82, Philosophy F). b2 = 0 cells are untouched (identical results).
  if (b2 == 1L) {
    aug <- .rec_augment_bspline(train, test, grid)
    train <- aug$train; test <- aug$test; grid <- aug$grid
  }

  model_plan <- .rec_model_plan(scenario_row)
  list(
    train = train, test = test, grid = grid,
    truth_grid = truth_grid,
    model_plan = model_plan,
    scenario_meta = list(
      scenario_id = scenario_row$scenario_id,
      regime = regime, bitmask = scenario_row$bitmask,
      pathologies_active = scenario_row$pathologies_active,
      gdpar_family = scenario_row$gdpar_family,
      K = scenario_row$K, rep = rep,
      truth_seed = truth_seed, realiz_seed = realiz_seed,
      n = n, n_train = n_train
    )
  )
}

## Simulate the response given the structural location predictor.
.rec_simulate_response <- function(regime, eta_loc, cov, b3, b4, b5) {
  C <- .rec_const
  n <- length(eta_loc)
  if (regime == "C") {
    sigma <- .rec_sigma_cont(cov, b4)
    if (b3 == 1L) {
      ## AR(1) errors with the target marginal SD (unit-variance series
      ## scaled per-observation by sigma).
      e <- .rec_ar1_series(n, C$ar_rho, innov_sd = 1)
      err <- e * sigma
    } else if (b5 == 1L) {
      ## heavy tails: scaled Student-t with unit-scale * sigma.
      err <- sigma * stats::rt(n, df = C$nu_heavy)
    } else {
      err <- stats::rnorm(n, sd = sigma)
    }
    if (b5 == 1L && b3 == 0L) {
      ## already heavy-tailed above
      return(eta_loc + err)
    }
    if (b5 == 1L && b3 == 1L) {
      ## heavy-tailed AR(1): AR series of t-innovations
      e_t <- numeric(n)
      e_t[1L] <- stats::rt(1L, df = C$nu_heavy)
      for (t in 2:n) e_t[t] <- C$ar_rho * e_t[t - 1L] +
        stats::rt(1L, df = C$nu_heavy)
      return(eta_loc + sigma * e_t)
    }
    return(eta_loc + err)
  }
  ## count regime
  lat <- if (b3 == 1L) .rec_ar1_series(n, C$ar_rho, innov_sd = 0.3)
         else rep(0, n)
  lambda <- exp(eta_loc + lat)
  if (b5 == 1L) {
    mu_counts <- stats::rnbinom(n, size = C$nb_size, mu = lambda)
  } else {
    mu_counts <- stats::rpois(n, lambda = lambda)
  }
  if (b4 == 1L) {
    pi <- .rec_pi_count(cov, b4)
    zero <- stats::rbinom(n, size = 1L, prob = pi)
    mu_counts[zero == 1L] <- 0L
  }
  as.numeric(mu_counts)
}

## True conditional moments E[y|x], SD[y|x], 0.95 quantile on the grid.
.rec_truth_moments <- function(regime, grid_cov, b2, b4, b5) {
  C <- .rec_const
  eta_loc <- .rec_eta_location(grid_cov, regime, b2)
  if (regime == "C") {
    sigma <- .rec_sigma_cont(grid_cov, b4)
    mu <- eta_loc
    if (b5 == 1L) {
      nu <- C$nu_heavy
      sd_y <- sigma * sqrt(nu / (nu - 2))
      tail_q <- mu + sigma * stats::qt(0.95, df = nu)
    } else {
      sd_y <- sigma
      tail_q <- mu + sigma * stats::qnorm(0.95)
    }
    return(list(mu = mu, sd = sd_y, tail_q095 = tail_q))
  }
  ## count regime: location is the poisson/nb rate lambda
  lambda <- exp(eta_loc)
  pi <- .rec_pi_count(grid_cov, b4)
  size <- if (b5 == 1L) C$nb_size else Inf
  base_mean <- lambda
  base_var  <- if (is.finite(size)) lambda + lambda^2 / size else lambda
  mu <- (1 - pi) * base_mean
  ## E[y^2] = (1-pi)*(var + mean^2); Var = E[y^2] - mu^2
  ey2 <- (1 - pi) * (base_var + base_mean^2)
  sd_y <- sqrt(pmax(ey2 - mu^2, 0))
  tail_q <- .rec_count_tail_quantile(lambda, pi, size, prob = 0.95)
  list(mu = mu, sd = sd_y, tail_q095 = tail_q)
}

## 0.95 quantile of a (zero-inflated) Poisson / Negative-Binomial, per
## grid point. Inverts the mixture CDF F(k) = pi*1{k>=0} + (1-pi)*F0(k).
.rec_count_tail_quantile <- function(lambda, pi, size, prob = 0.95) {
  n <- length(lambda)
  out <- numeric(n)
  for (i in seq_len(n)) {
    lam_i <- lambda[i]; pi_i <- pi[i]
    cdf0 <- function(k) {
      if (is.finite(size)) stats::pnbinom(k, size = size, mu = lam_i)
      else stats::ppois(k, lambda = lam_i)
    }
    k <- 0L
    repeat {
      f <- pi_i + (1 - pi_i) * cdf0(k)
      if (f >= prob || k > 1e6) break
      k <- k + 1L
    }
    out[i] <- k
  }
  out
}

## ---------------------------------------------------------------------
## Model plan: per-method modeling recipe derived from the scenario.
## gdpar uses the (formula, amm) path for K=1 and the gdpar_bf path for
## K>=2. Slot specifications are DGP-aware (B9.12): each slot models the
## covariate dependence the ground truth actually has, mirroring the
## logic the gaussian regime already uses (K=2 only when b4=1):
##   - location mu: x1,x2 (+ W(x1) if b2).
##   - continuous scale sigma: a(z) iff b4=1 (heteroscedasticity), else
##     intercept-only (truth has constant sigma when b4=0).
##   - student_t nu: intercept-only always (the factorial never modulates
##     the tail index with covariates; truth nu = 3 constant).
##   - count overdispersion phi (neg_binomial_2/zinb): intercept-only
##     always (truth nb_size constant).
##   - zero-inflation pi (zip/zinb): a(z) (truth pi = inv_logit(zi0+zi1*z)).
## Specifying a slot on a noise covariate (the old nu~a(w)/pi~a(w)) or on
## z where the truth is constant (the old phi~a(z)) is spurious
## over-parametrization: it broke convergence (student_t) and prevented
## recovery of the true structure (zinb pi modelled on noise w).
## ---------------------------------------------------------------------

.rec_model_plan <- function(scenario_row) {
  regime <- scenario_row$regime
  b2 <- scenario_row$b2
  fam <- scenario_row$gdpar_family
  K <- scenario_row$K
  list(
    regime = regime,
    gdpar_family = fam,
    K = K,
    use_W = (b2 == 1L),                 ## nonlinearity stratum active
    ## D82 / Philosophy F (refined by D83): native smooth of x1 per method.
    ## mgcv/brms use s(x1); INLA uses f(x1, model="rw2") (D83); gdpar uses
    ## a(bs(x1)); rstanarm (no native smoother) consumes these frozen shared
    ## basis columns. gdpar and rstanarm both read smooth_cols; INLA ignores
    ## them and bins x1 itself (inla.group) for its rw2 latent model.
    smooth_cols = if (b2 == 1L) .rec_smooth_cols() else NULL,
    mean_vars = c("x1", "x2"),
    scale_var = "z",                    ## 2nd slot driver (sigma / phi / pi)
    scale_hetero = (scenario_row$b4 == 1L), ## continuous scale depends on z
    response_type = if (regime == "C") "continuous" else "count"
  )
}

## Convenience: the full 64-scenario list keyed by scenario_id, each a
## closure realizing rep. Mirrors scenarios_synthetic_hard_multi() shape
## for driver familiarity, but parameterized by rep.
scenarios_recovery_factorial <- function() {
  tab <- recovery_scenario_table()
  out <- lapply(seq_len(nrow(tab)), function(i) {
    row <- tab[i, , drop = FALSE]
    function(rep = 1L, ...) recovery_realize_cell(row, rep, ...)
  })
  names(out) <- tab$scenario_id
  out
}

# gdpar 0.1.0 (first public release, 2026-06-27)

First public release of `gdpar`. It implements the General Dynamic Parameter
framework, in which each unit's parameter is decomposed as
`theta_i = theta_ref + Delta(x_i, theta_ref)` around a population reference, with
`Delta` following the canonical Additive--Multiplicative--Modulated (AMM) form
`a(x) + b(x) * theta + W(theta) x`. Estimation is via Path 1 (hierarchical
Bayesian inference through 'Stan' / 'cmdstanr'); Paths 2 (varying-coefficient)
and 3 (hypernetwork) are specified at reference grade and queued for future
versions.

Packaging note (pre-CRAN correction): the fitting entry points now run all
input validation before the optional `cmdstanr` dependency is required, which
is checked only immediately before it is used. The package and its full test
suite therefore work on machines where the `Suggests` package `cmdstanr` is not
installed. This is a packaging-only change with no effect on results when
`cmdstanr` is available.

This release consolidates the development cycles documented in detail below
(Sessions 8.3--8.6, the Block 8.4 audit, and the Block 9 re-validation):

* **Distributional regression, `K >= 1` slots.** Gaussian, Student-t, Gamma,
  Beta, Tweedie, Poisson, negative-binomial, zero-inflated and hurdle families,
  plus heterogeneous per-slot families, with B-spline `W` bases and arbitrary
  covariate dimension `p`.
* **Empirical-Bayes and full-Bayes estimation**, with EB-vs-FB comparators
  (marginal total variation and joint KSD) and a unified Path C bit-exact
  comparator for the `K > 1` and `p > 1` regime.
* **Causal positioning.** A T-learner CATE/ITE bridge on the AMM side and an
  external meta-learner comparator (grf, EconML via reticulate).
* **Geometric-robustness sampling engine (Block RG).** An opt-in,
  default-bit-identical Riemannian / sub-Riemannian / Finsler metric stack with a
  certifying orchestrator and a learned expected-Fisher metric, for stiff and
  near-deterministic posteriors.
* **Dependence-robust inference (Block 9, Axis 2).** Diagnostics (lag-1 /
  Durbin-Watson; Moran's I) and block-bootstrap robust standard errors for
  temporally and spatially dependent data, with a data-driven (Politis-White)
  block length, on both the Empirical-Bayes and full-Bayes paths. `gdpar` does
  not model the dependence; only its inference is made robust to it.
* **Identifiability diagnostics and validity tests** for the population
  reference, residual diagnostics (Dunn-Smyth quantile residuals; G1/G2/G3), and
  a golden-regression test suite.

**Validation.** Block 9 closed adversarial internal re-validation, a synthetic
recovery benchmark against mgcv / brms / INLA / rstanarm with known ground truth
(3200 cells), and an organic eBird benchmark (80 cells). `gdpar` leads on the
distributional structure it models (heavy-tail quantiles, zero-inflated count
means, heteroscedastic scale), is statistically indistinguishable from the
strongest competitor (mgcv) on generic predictive accuracy, is robust where it
makes no modelling claim (e.g. autocorrelation), and is the most computationally
expensive method in the roster. Reports under `inst/benchmarks/results/`.

Detailed per-session development notes follow.

## Axis 2 — full-Bayes dependence-robust inference (Session B9.39, 2026-06-26; D102)

* **The dependence diagnostics and the block-bootstrap robust standard errors now
  accept a scalar full-Bayes fit (`gdpar_fit`), not only a scalar Empirical-Bayes
  fit (`gdpar_eb_fit`).** `gdpar_dependence_diagnostic()`,
  `gdpar_dependence_robust()`, `gdpar_spatial_dependence_diagnostic()` and
  `gdpar_spatial_dependence_robust()` are unchanged in name and signature; they
  simply work on either path now, closing the EB/FB asymmetry left open at the
  start of Axis 2. The K > 1 / p > 1 paths remain deferred on both EB and FB.
* **One shared engine, three class-dispatched touchpoints.** The refit loop, the
  Politis-White / spatial block-length selectors, the temporal and spatial
  diagnostics and the honesty contract are path-agnostic; only the per-fit
  extraction of the point estimate, the model SE and the Dunn-Smyth residuals is
  dispatched by object class. The Empirical-Bayes code paths are byte-identical,
  so the existing bit-exact EB regression gate is preserved.
* **Full-Bayes point estimate / model SE = posterior mean / posterior SD** of
  each AMM coefficient (`theta_ref`, `a_coef`, `b_coef`, `W_raw`, the last on its
  raw scale, for parity with the Empirical-Bayes extractor; the `theta_ref`
  hyperparameters are excluded). `robust_se` is the block-bootstrap SD of the
  per-refit posterior means and `se_ratio = robust_se / model_se` remains a
  like-for-like (SD-vs-SD) ratio. Each refit re-runs the full HMC and is
  markedly more expensive than an Empirical-Bayes refit (the opt-in cost message
  says so).
* **Per-refit convergence accounting** (new `refit_diagnostics` field on the
  robust result: `max_rhat`, `min_ess_bulk`, `n_divergent_refits`,
  `n_high_rhat_refits`). Refits are never excluded or down-weighted (that would
  be a non-random, biasing screen); a single informational note is emitted only
  when a refit's R-hat clearly exceeds 1.05. This applies to both paths.
* **Documented full-Bayes caveats:** the finite-iteration Monte-Carlo error of
  each refit's posterior mean conservatively (slightly) inflates `robust_se`;
  under an informative prior `se_ratio < 1` is benign prior regularization, not
  an overstated model SE; and the Empirical-Bayes `theta_ref` point estimate (the
  Laplace mode) and the full-Bayes one (the posterior mean) are different
  estimands that coincide only asymptotically. A bagged / widened posterior
  (BayesBag) is recorded as a deferred lateral, not adopted. The full-Bayes
  design was cross-lineage-reviewed; the review's concrete proposals (a hard
  minimum-ESS abort, an estimator-replacing bag) were audited and **not** adopted
  in favour of the more robust informational accounting, while its documentation
  findings (the `se_ratio < 1` and mode-vs-mean subtleties) were adopted.

## Axis 2 — data-driven block length (Session B9.38, 2026-06-26; D101)

* **`gdpar_dependence_robust()` and `gdpar_spatial_dependence_robust()` gain an
  opt-in data-driven block size, `block_length = "auto"` / `block_size =
  "auto"`.** The default behaviour is unchanged (the rate-only `n^(1/3)` /
  `n^(1/4)`); `"auto"` replaces the rate's *constant* with a value selected from
  the fitted residuals, with no extra Stan refit, and the rate as the fallback.
  Both functions also gain `residual_type` / `randomize_seed`, used only to feed
  the `"auto"` selector. The chosen value and method are reported
  (`block_length` / `block_size` and the new `block_length_method` /
  `block_size_method`).
* **Temporal: the canonical Politis & White (2004) selector** (with the Patton,
  Politis & White 2009 correction), hand-rolled in base R (no **np**
  dependency): the adaptive flat-top spectral estimate of the residual
  autocovariances gives `b_opt = (2 ghat^2 / D)^(1/3) n^(1/3)` with the
  overlapping moving/circular-block constant `D = (4/3) spec^2`; longer
  dependence yields a longer block, white noise yields unit blocks, and a
  degenerate series falls back to the rate. A `p << n` caveat is documented.
* **Spatial: a data-driven calibration over the cell grid** (Politis & White has
  no established spatial plug-in). For each candidate `g`, cheap (no-refit)
  spatial block resamples give the bootstrap variance of the design-weighted
  residual functionals (the coefficient's influence directions); `g` minimises
  an empirical mean-squared error — squared bias anchored at the *largest* blocks
  (least biased) plus a `1/n_tiles` sampling-variance term (Lahiri 2003). A
  decorrelating cross-lineage review supplied the MSE skeleton; two of its
  choices were corrected after an audit and empirical validation (the bias
  anchor, and the variance term), because the proposed forms would have made the
  selector either anticonservative or non-adaptive. A single isotropic `g` is
  used (anisotropy is a documented, deferred limitation).

## Axis 2 — spatial dependence-robust inference (Session B9.37, 2026-06-26; D100)

* **`gdpar_spatial_dependence_diagnostic()` (new, exported) — Moran's I for the
  residuals of a scalar Empirical-Bayes fit.** The spatial sibling of
  `gdpar_dependence_diagnostic()`: it builds a row-standardized spatial weight
  matrix (default `k`-nearest-neighbour, the data-driven `k = max(4, min(round(
  log n), n - 1))`; or a distance band; or a user-supplied `W` carrying domain
  knowledge), hand-rolls Moran's I in base R (no **spdep** / **sf** dependency),
  and reports a significance test — a two-sided permutation test by default
  (robust to non-normal Dunn-Smyth residuals and asymmetric weights), or the
  cheaper analytic Cliff-Ord normal approximation (which warns under an
  asymmetric `W`). Guards cover isolated locations (zero-weight rows warn and
  return `NA`), duplicate coordinates (permitted), small `n` (hard/soft
  warnings), and the lon/lat (project-first) and model-misspecification caveats
  are documented. gdpar does **not** model the spatial dependence; this only
  makes the violation visible.

* **`gdpar_spatial_dependence_robust()` (new, exported) — spatial
  block-bootstrap-by-refit standard errors.** Refits the model on `B` spatial
  block-bootstrap resamples and reports bootstrap SEs and percentile intervals
  alongside the model-based ones, in the working-independence + robust-variance
  spirit of Liang & Zeger (1986) (point estimates unchanged). Default scheme is
  non-overlapping **tiled** blocks with a **randomized grid origin** per
  replicate (Politis-Romano-Lahiri), plus an opt-in overlapping **moving**
  scheme. The default block side is `g = max(2, round(n^(1/4)))` cells per axis,
  the `d = 2` case of the variance-MSE-optimal rate `M ~ n^(d/(d+2))` points per
  block, which reduces **exactly** to the temporal `n^(1/3)` block length at
  `d = 1` (derivation and the registered cross-lineage dissent toward `n^(1/6)`
  are in `?gdpar_spatial_dependence_robust`). Collinear coordinates abort; a
  single-cell collapse warns.

* **Internal refactor (regression-gated).** The temporal `gdpar_dependence_
  robust()` and the new spatial function share one extracted engine
  (`.gdpar_dependence_robust_engine`) and one residual extractor
  (`.gdpar_dependence_residuals`); the temporal path is verified **bit-identical**
  pre/post refactor by a fixed-seed regression gate (a frozen pre-refactor copy
  vs. the live engine-backed function, same robust-SE table) plus a Stan-free
  golden that locks the engine's RNG-consumption order. No change to the
  temporal API or output.

## Bloque RG — RG.7 (step 4): the automated Laplace fallback (Session B9.35, 2026-06-25; D99)

* **`gdpar_geom_laplace()` (new, exported) — the Laplace approximation as a
  first-class capability.** Given a `gdpar_geom_target`, it climbs to the
  posterior mode (L-BFGS-B warm start + a modified-Newton polish when an exact
  Hessian is exposed, reading the **same** target and gradient the sampler uses),
  forms the precision `M = -`Hessian (the observed information, exact cmdstan
  Hessian or finite differences), reports its covariance `M^{-1}`, optional iid
  draws `mode + M^{-1/2} z`, and — crucially — a **fidelity diagnostic** of the
  Gaussian against the true posterior over the same draws: the
  importance-sampling ESS, the PSIS Pareto-`k` (when **loo** is installed), and
  the mean/max log-density drop against its Gaussian expectation `d/2`. These are
  distilled into a single scalar label `"good"` / `"poor"` / `"very_poor"` so the
  approximation is never mistaken for exact MCMC. The **un-floored** condition
  number is reported and a saddle (non-positive-definite curvature) is flagged
  loudly. This promotes the validated extractor core of RG.7
  (`inst/benchmarks/scripts/rg7_laplace_elpd.R`) into `R/`.

* **`gdpar_geom_orchestrate()` and `gdpar_geom_fit()` gain `laplace_fallback =
  FALSE` (opt-in) and `laplace_draws = 0L`.** When opted in and a run ends in a
  **certified limit** (the genuinely non-Gaussian canyon the geometry ladder
  cannot sample — RG.7's eBird Tweedie count), the orchestrator attaches the
  Laplace approximation (`$laplace`) and relabels the status
  `"certified_limit_laplace"`: the sampling limit still stands **and** a labelled,
  fidelity-diagnosed Laplace posterior is provided — exactly the mgcv/REML and
  INLA/Laplace competitor-parity regime. On the out-of-scope path the fallback is
  attached only when the curvature at the mode is genuinely positive-definite
  (status `"out_of_scope_laplace"`); otherwise the certificate is left untouched
  (no overreach). **The default `FALSE` leaves the output bit-identical** (no
  `$laplace`, the original status set), so the goldens stay bit-identical.

## Bloque RG — RG.6 (part ii): integration of the geometry-adaptive orchestrator + the pedagogical vignette (Session B9.31, 2026-06-03)

* **`gdpar_geom_bridge()` (new, exported) — the durable, path-agnostic core.**
  Turns an already-fitted `gdpar` object into the inputs the geometry-adaptive
  controller `gdpar_geom_orchestrate()` consumes, **without touching `gdpar()`**.
  It reads the fit's compiled `cmdstan` model and Stan data, exposes the
  standalone `log_prob` / `grad_log_prob` / `hessian` methods, derives the
  unconstrained dimension and a posterior-mean warm-start, and returns a
  `(target, geom_target, fisher, reference)` tuple. Because the diagnostic needs
  a *re-samplable* model and a `CmdStanMCMC` object cannot be re-sampled, the
  bridge recompiles one from the fit's own Stan source (`$code()`, with the
  standalone methods; cmdstanr's content-hash cache makes it a cache hit), while
  the engine target reuses the fitted object directly. The Hessian compilation is
  **best-effort**: if higher-order autodiff will not compile (models with custom
  densities, the Tweedie `lpdf` being the case RG.7 targets), it falls back to
  gradient-only methods, which still serve the Euclidean, dense and
  sub-Riemannian (expected-Fisher) levels. This is the tool RG.7 points at the
  real Tweedie count of benchmark 9.2.O.

* **`gdpar_geom_fit()` (new, exported) — the one-call ergonomic entry.** A
  *sister* of `gdpar()` (not an internal branch) for the K-individual path: it
  builds and compiles the model through the **shared seam `.gdpar_K_build()`**,
  then runs the orchestrator instead of the default NUTS fit, returning a rich
  `gdpar_geom_fit` object. Print methods for both new classes are provided.

* **`.gdpar_K_build()` (internal) extracted from `.gdpar_K()`.** A
  behaviour-preserving refactor of the model-building phase (everything up to but
  not including `cs_model$sample`): `.gdpar_K()` and `gdpar_geom_fit()` now share
  a **single source** (no duplication, no throwaway computation). The default
  branch (`compile_model_methods = FALSE`) issues the byte-identical
  `cmdstanr::cmdstan_model(stan_path)` call used before, so the K-path goldens
  stay bit-identical (verified: the K=2 Gaussian golden re-fits with
  `max|delta| = 0`); `gdpar_geom_fit()` passes `compile_model_methods = TRUE` to
  expose the standalone methods.

* **Pedagogical vignette `vop08_geometric_robustness.Rmd` (new).** An abundant,
  two-level narrative of the whole Block RG (RG.0 → RG.6): the conceptual bridge
  to the user's geometry manuscript (with an organic critical reading of its §12
  hierarchy and its §14 "future statistical line", and of the ORPHEUS §16.3
  demonstrated-vs-conjectured honesty), the diagnostic with size-invariant
  signals, the five sampler levels, the orchestrator with its armour and its
  certified limit, Option A, and the integration. The **evaluated** chunks are
  cheap (suite metadata, the budget/criteria/thresholds, and the geometry engine
  on pure-R closure targets — no Stan); the heavy `cmdstan`-backed runs are shown
  with `eval = FALSE` and reproduced end to end by the gated script
  `inst/scripts/geometry_pilots_deep.R` (`GDPAR_RUN_GEOMETRY_PILOTS=1`).

  Opt-in throughout: the default `gdpar()` fit path is untouched, so the goldens
  remain bit-identical (canonized as **D97**, session B9.31).

## Bloque RG — RG.6 (part 1): Option A folded in — `sigma_a_k` compaction (Session B9.30, 2026-06-03)

* **Compaction of the per-slot additive scale `sigma_a_k` to the slots that
  carry free `a` coefficients (`J_a_free > 0`).** A distributional / multivariate
  slot whose `a()` component contributes no free coefficient — either because the
  slot declares no `a()` at all (`use_a_k == 0`, e.g. `phi ~ 1`, `p ~ 1`) or
  because it declares an intercept-only `a()` (`a = ~ 1`, fully absorbed into the
  anchor `theta_ref`) — used to declare a sampled-but-unused half-prior scale
  `sigma_a_k[k]`. That scale is a **flat direction** (a non-identified nuisance
  dimension) that the no-U-turn sampler explores pointlessly and that, combined
  with other pathologies, produced the `rhat = Inf` first diagnosed in session
  B9.21. The Stan templates now compute, in `transformed data`, the count
  `n_sigma_a` of slots carrying free `a` coefficients and an index map
  `sigma_a_idx` (the compacted position per slot, `0` if none), declare
  `sigma_a_k` with length `n_sigma_a`, and index it as `sigma_a_k[sigma_a_idx[k]]`
  wherever the additive coefficients are reconstructed or given their hierarchical
  prior. Touched: the three K-individual canonical pieces
  (`amm_canonical_distrib_K.stan`, `amm_canonical_eb_conditional_K.stan`,
  `amm_canonical_eb_marginal_K.stan`) via the shared codegen `generate_a_blocks_K`,
  the three K×p multivariate pieces (`amm_canonical_pmulti_KxP.stan`,
  `amm_eb_conditional_KxP.stan`, `amm_eb_marginal_KxP.stan`) inline, and the Path C
  K×p random-init helper so the init vector matches the compacted length.

  Methodological decision (session B9.30, the cornerstone rule of maximum
  multidimensional robustness; canonized as **D96**):

  - **The criterion is `J_a_free > 0`, not `use_a_k == 1`.** The mathematically
    correct definition of when `sigma_a_k` exists is "there is at least one free
    `a` coefficient for it to scale". This single criterion removes the flat
    direction in **both** sub-cases (no `a()` and intercept-only `a()`), so it is
    strictly more robust than a `use_a_k`-based rule, which would miss the
    intercept-only case. Normalising the inconsistency upstream in `amm_spec`
    (forcing `use_a_k = 0` for an intercept-only `a()`) was rejected as an
    over-correction: "an `a()` was declared" and "an `a()` has free coefficients"
    are legitimately distinct facts, and the truth that governs `sigma_a_k` lives
    at the Stan level, where `J_a_free > 0` states it exactly.

  - **Bit-identity and the one re-bootstrapped golden.** When every slot carries
    free `a` coefficients (the case of every golden whose formula puts `a(x)` on
    each slot) `n_sigma_a == K` and `sigma_a_idx` is the identity, so the generated
    model is mathematically identical and the frozen golden draws stay bit-exact
    (verified by re-fit-and-compare on the K, K×p and EB paths). The sole
    exception is the K×p `student_t` golden `g_KxP_004`, whose `nu` slot declares
    an intercept-only `a()` (`a = ~ 1`) and therefore carried a flat
    `sigma_a_k[3]`; Option A removes it. That golden was re-bootstrapped, and the
    change was validated as a **faithful marginalisation**: `sigma_a_k[3]` appears
    only in its own prior (never in the likelihood or any other parameter's
    prior), so it is a-priori independent of the rest and its removal leaves the
    joint over the surviving parameters mathematically invariant; the re-fit shows
    exactly one fewer parameter and no surviving marginal mean disagreeing beyond
    Monte-Carlo error. A regression guard test fixes that a slot with an
    intercept-only `a()` yields a `sigma_a_k` of length `< K`.

## Bloque RG — RG.5: geometry-adaptive sampling orchestrator (Session B9.29, 2026-06-03)

* **Geometry-adaptive orchestrator (`gdpar_geom_orchestrate()`).** The closed loop
  of the Block RG charter: an opt-in controller that diagnoses the geometry of a
  posterior (reusing `gdpar_geometry_diagnostic()`), classifies the pathology,
  selects a level of the sampler hierarchy (`euclidean_diagonal` →
  `euclidean_dense` → `riemannian` → `relativistic` → `sub_riemannian`), samples
  with that level (reusing `gdpar_geom_hmc()` and the RG.2/RG.3/RG.4 metrics),
  re-diagnoses, and either escalates the level or emits a **certified limit**. It
  is standalone and does not touch the package's fit path, so the default branch
  is bit-identical (goldens intact). Returns an object of class
  `gdpar_geom_orchestration` (`status` one of `"resolved"`, `"certified_limit"`,
  `"out_of_scope"`), carrying the per-round decision/sampling ledger, the winning
  metric, the best result so far, the budget spent and a reproducibility block.

  Methodological decisions of this sub-phase (session B9.29, four hinges consulted
  upfront and decided by the user under the cornerstone rule of maximum
  multidimensional robustness):

  - **Level selection — the combined map.** The diagnostic's transparent,
    calibrated rule-based classifier is the primary selector; a continuous
    proximity score over the size-invariant signals is a second, independent
    estimator. When they agree and confidence is high, the discrete level is used
    directly; when they conflict or confidence is low, the controller starts at
    the lower of the two candidate levels and lets the escalation climb — robust
    at the funnel/heavy-tail border the RG.1.c calibration found mutually
    confusable. A user `level_map` / `entry_level` overrides this.

  - **Escalation — the closed-loop, re-diagnosis-guided controller, fully
    armoured.** On a failed level the controller re-diagnoses with a fresh pilot
    and lets the new signature pick the next level (a re-diagnosis-guided upward
    jump that may skip ladder rungs),
    under a monotone ratchet (never below a level already tried), a visited-state
    memo (a provably acyclic walk), a global round cap, per-level caps, a
    no-progress (stall) detector, hierarchical budget accounting with cost-aware
    admission and cheap-probe-then-full tiering (the successive-halving spirit; Li
    et al. 2018), a coarse step-size search per level (so a level is not failed for
    a mis-set step), a per-fit wall-time watchdog, graceful degradation (the best
    result so far is always returned), deterministic seeding (the adaptive
    trajectory is bit-reproducible, so a re-run retraces it — the resumability
    guarantee, with an optional atomic ledger checkpoint), a multi-signal success
    gate with hysteresis, a frozen-level final sampling phase (preserving the
    exactness / ergodicity of the returned draws; cf. the diminishing-adaptation
    condition of Roberts & Rosenthal 2007), and the capability-subsumption
    guarantee (level L+1 strictly generalises level L, so the monotone ladder
    fallback is a proven monotone-improvement backstop when an adaptive jump
    misfires). Every level is Metropolis-exact, so a wrong jump wastes only a
    bounded amount of budget, never correctness.

  - **The certified limit (`gdpar_geom_certificate`).** When the budget is
    exhausted without success the controller returns a first-class, reproducible,
    falsifiable certificate: the **demonstrated** evidence in three rigour layers
    (algebraic = the geometry; statistical = the per-level sampler diagnostics;
    numerical = the budget and fits) plus a **conjectured prescription** — the
    smallest change predicted to break the limit, tagged as a conjecture in the
    sense of the demonstrated/conjectured honesty convention and made falsifiable
    by a re-run test. A diagnosis pointing outside the geometry ladder
    (multimodality, a boundary, a flat direction) short-circuits to a certificate
    naming the proper remedy (tempering, reparametrisation, Option A), without
    overreaching by pretending a geometry metric fixes it.

  - **Scope — standalone with a rich return.** The controller never touches
    `gdpar()`; the return carries everything a future integration (RG.6) will need
    (the winning metric, the ledger, the reproducibility block) without freezing a
    contract that the integration might have to redesign. The statistical-layer
    validation exercises a real compiled cmdstan model end to end.

  The tunable factories `gdpar_geom_orchestrate_budget()` (budget and stopping
  rule) and `gdpar_geom_orchestrate_criteria()` (the multi-signal success gate)
  expose the controller's data so it can be re-calibrated. Anchors: Betancourt
  (E-BFMI, sampling pathologies); Girolami & Calderhead 2011; Roberts & Rosenthal
  2007 (diminishing adaptation / containment); Li et al. 2018 (Hyperband /
  successive halving); Xu et al. 2008 (SATzilla, runtime-predicted algorithm
  selection with a presolver fallback). gdpar integrates these established
  patterns; it does not claim to invent them.

## Bloque RG — RG.4 (cont.): Finsler / relativistic geometry, the heavy-tail remedy (Session B9.28, 2026-06-03)

* **Finsler / relativistic metric (`gdpar_geom_metric_relativistic()`).**
  The level-4 geometry of the Block RG hierarchy and the remedy for heavy tails
  and directional anisotropy (the `G3_heavy_tails` target, *not* the count). A
  Gaussian kinetic energy gives the unbounded velocity `M^{-1} p`, so a large
  momentum drawn in a heavy tail produces an arbitrarily large position step
  that a fixed-step integrator cannot follow (overshoot, divergences, a step
  size forced to the stiffest region). The relativistic kinetic energy caps the
  velocity at a finite speed `c`, taming the tails and the ill-conditioning
  while staying exact. Coupled to the position-dependent Riemannian metric
  `M(theta)` of RG.3 (the chosen, maximally robust form), the kinetic energy is
  the relativistic energy of a particle of rest mass `m` on the statistical
  manifold,
  `K = c sqrt(p' M(theta)^{-1} p + m^2 c^2) + 0.5 log det M(theta)`.

  Methodological decisions of this sub-phase (Lu, Perrone, Hasenclever, Teh and
  Vollmer 2017 for the relativistic kinetic; Livingstone, Faulkner and Roberts
  2019 for the kinetic-energy choice under heavy tails; Girolami and Calderhead
  2011 for the generalised implicit leapfrog; Randers 1941 for the asymmetric
  metric), with an organic-critical reading of the user's geometry document
  (section 12.3):

  - *Bounded velocity by construction.* The velocity
    `grad_p K = c M^{-1} p / sqrt(p' M^{-1} p + m^2 c^2)` has `M`-norm strictly
    below `c` for every momentum — the property that domps the tails. The
    `0.5 log det M` term is the same normaliser as the Gaussian Riemannian
    kinetic, so the `theta`-marginal of the joint is exactly the posterior: the
    kinetic energy is a preconditioner, not part of the target, and the
    Metropolis correction with the exact density keeps the sampler exact for any
    `speed` and `rest_mass` (they govern only efficiency).

  - *Strict generalisation of the Riemannian level.* As `c -> infinity` the
    kinetic energy reduces to the Gaussian Riemannian kinetic of RG.3 (the
    non-relativistic limit); larger `speed` mixes faster in the bulk with less
    tail-taming, smaller `speed` caps sooner and is more robust.

  - *Finsler structure, with the asymmetric Randers piece deferred on purpose.*
    A relativistic kinetic energy is the Legendre dual of a Finsler norm on
    velocities — a norm not induced by an inner product, the speed `c` being the
    Finsler unit ball. The asymmetric Randers extension `F = sqrt(g(v,v)) +
    beta(v)` (an irreversible drift) is **not** included: it makes the kinetic
    energy odd in `p` and would break the reversibility the Metropolis
    correction relies on. It models irreversible *dynamics* — the document's
    ontological interest, the cost of raising complexity differing from lowering
    it — a different goal than exact *sampling* of a fixed posterior. The even,
    bounded relativistic kinetic realises exactly the part of the Finsler
    insight that serves robust, exact sampling.

  - *Non-separable Hamiltonian -> dedicated integrator.* Because `K` depends on
    both `theta` (through `M`) and `p`, the Hamiltonian is non-separable: the
    explicit leapfrog does not apply, and the existing implicit leapfrog
    hardwires the Gaussian velocity. A dedicated generalised implicit leapfrog
    (three reversible, volume-preserving sub-steps with the relativistic
    velocity) is carried in the metric's `integrator` slot, run by
    `gdpar_geom_hmc()` in place of the default leapfrog; the default branch stays
    bit-identical. The momentum is refreshed from the exact relativistic
    momentum law (an inverse-CDF radial sampler under `p = L q`).

  Opt-in throughout: the default branch is bit-identical, no fit path is touched
  and the golden regression is intact. `gdpar_geom_hmc()` now uses a metric's
  own kinetic energy when it carries one (and rejects, rather than crashes on, a
  proposal that drives the kinetic or target undefined). Validated in three
  layers — algebraic (kinetic gradients vs finite differences, the `M`-norm
  velocity bound, the Gaussian limit, the dedicated integrator's reversibility
  and energy conservation, the momentum law); statistical ungated (a heavy-tailed
  Student-t recovered); statistical gated (a correlated multivariate-t where the
  Euclidean kinetic overshoots, and a real cmdstan-backed correlated
  multivariate-t with a position-dependent SoftAbs mass matching the no-U-turn
  sampler). Decision **D94**.

## Bloque RG — RG.4: sub-Riemannian geometry, the quasi-deterministic remedy (Session B9.27, 2026-06-03)

* **Sub-Riemannian metric and integrator (`gdpar_geom_metric_subriemannian()`).**
  The level-5 geometry of the Block RG hierarchy and the remedy for a
  quasi-deterministic posterior (the eBird count / tweedie case), where the
  typical set contracts onto a lower-dimensional manifold: the expected Fisher
  information grows without bound along the stiff "wall" directions while the
  "floor" directions still carry genuine variation. A standard or even
  Riemannian sampler is then forced to a vanishing step by the walls; this
  geometry equips only a distribution of accessible directions with motion and
  glides along the floor instead.

  Methodological decisions of this sub-phase (Montgomery 2002; Shahbaba,
  Lan, Johnson and Neal 2014 for the split):

  - *Accessible distribution from the near-null space of the expected Fisher.*
    The Fisher (closed-form, or estimated by `gdpar_geom_fisher_simulator()`) is
    eigendecomposed at a reference; a **continuous** verticality filter
    `w_i = sigma((log lambda_i - log tau) / softness)` blends each direction
    between floor (small `lambda`) and wall (large `lambda`) with no hard cut, so
    a borderline direction is split smoothly and reversibility is never broken by
    a discontinuous classification.
  - *Exact wall flow, no step-size penalty.* The wall curvature
    `A = U diag(w_i lambda_i) U^T` defines a fixed reference quadratic whose
    harmonic flow is integrated in closed form (a symplectic rotation per mode);
    the stiff walls cost no step size while the floor follows the free-drift
    limit. A Strang splitting (half residual kick, exact flow, half kick) keeps
    the scheme symplectic and time-reversible. The default threshold `tau` is the
    floor scale, which caps the leapfrog residual at the floor so the step is
    governed by the floor and never the walls; `suggested_epsilon` exposes that
    step. The Metropolis correction with the exact density keeps the sampler
    exact however coarse the Gaussian treatment of the walls is — only
    efficiency, never correctness, depends on it.

  Validated in three layers: algebraic (the exact flow conserves the reference
  quadratic and reverses to machine precision; the Strang trajectory is
  reversible; exact on a pure quadratic at any step; the filter is monotone and
  caps the residual); statistical (a mild canyon recovered through
  `gdpar_geom_hmc`); and gated heavier runs — the G4 canyon at `n = 1000` (the
  sub-Riemannian sampler accepts 0.90 with zero divergences and recovers both
  scales where the Euclidean sampler diverges on every proposal at the same
  step) and a **real cmdstan-backed near-deterministic Poisson count** (the
  simulation-based Fisher matches the analytic information to 2%, and the
  sampler matches NUTS in mean and standard deviation at a step that wrecks the
  Euclidean sampler) — the structural rehearsal for the eBird tweedie of RG.7.

* **Hardened proposal evaluation in `gdpar_geom_hmc()`.** A proposal that drives
  the target undefined (a cmdstan log-density that throws on a non-finite
  unconstrained value, an overflowing gradient) is now caught and counted as
  divergent rather than crashing the run. Stable proposals never trigger the
  handler, so the default branch stays bit-identical.

## Bloque RG — RG.3 Capa 2 completed: simulation-based Fisher + active learning (Session B9.26, 2026-06-03)

* **Simulation-based estimator of the expected Fisher information
  (`gdpar_geom_fisher_simulator()`).** Completes the learned Riemannian metric
  for models with no closed-form Fisher: the expected Fisher
  `I(theta) = E_y[ s s^T ]` is estimated by the average outer product of the
  log-likelihood score over data sets simulated from the model at `theta`,
  positive semi-definite by construction and unbiased. The estimate is a
  deterministic, reproducible function of `theta` (the RNG is reseeded from a
  position key and restored). The SoftAbs mean of the Gaussian-process surrogate
  acts as a structural control variate, so few replicates per site suffice. It
  plugs into the `fisher` slot of `gdpar_geom_metric_gp_fisher()`. Validated
  against closed forms: a constant Fisher (a normal location model) and a
  position-dependent one (a Poisson log-rate model), to a few percent.

* **Generative targets.** `gdpar_geom_target()` gains optional `simulate(theta)`
  and `score(theta, y)` slots (the per-data-set log-likelihood score whose outer
  product is the Fisher), distinct from the fixed-data `grad_log_prob`. Default
  branch unchanged.

* **Adaptive Riemannian HMC with online active learning
  (`gdpar_geom_rmhmc_adaptive()`).** The full novelty-driven loop: the learned
  metric is held fixed within each trajectory (preserving reversibility) and
  re-learned only between rounds from a reservoir that grows where the
  surrogate's epistemic novelty is high, with a decreasing (Robbins--Monro-style)
  number of admitted sites per round so the metric sequence settles; a final
  sampling phase uses the frozen metric. The sampler is **exact** in every phase
  regardless of the metric (the metric is a preconditioner, not part of the
  target), so no delayed acceptance is needed; only the efficiency is heuristic
  and is measured (E-BFMI, acceptance, novelty trace), not asserted (ORPHEUS
  section 16.3). A gated end-to-end demonstration runs the full pipeline over a
  real cmdstan-backed Poisson count model — the structural rehearsal for the
  eBird Tweedie of RG.7.

## Bloque RG — RG.3 Capa 2: learned expected-Fisher metric (Session B9.25, 2026-06-02)

* **A learned, smooth surrogate of the expected Fisher information
  (`gdpar_geom_metric_gp_fisher()`).** The general realisation of the natural
  Rao–Amari Riemannian metric where no closed form exists, building
  `M(theta) = L(theta) L(theta)^T` from a Gaussian-process surrogate of the
  log-Cholesky factor. Following an organic-critical reading of ORPHEUS-PIMC-A
  (the new sections 8 and 16), the surrogate's **mean function is the SoftAbs
  curvature** of Capa 1 and the process learns only the **smooth residual** to
  the expected Fisher at the reservoir sites. A single object unifies three
  ORPHEUS components without their drawbacks: (i) the metric degrades
  **continuously** to the always-available SoftAbs far from the reservoir (the
  kernel decays, the posterior mean returns to its SoftAbs mean) — no hard
  metric switch that would break reversibility; (ii) the predictive variance is
  the principled **epistemic-novelty** detector (`novelty(theta)`); (iii)
  positive-definiteness, smoothness and an analytically differentiable
  `dmass(theta)` hold by construction. The sampler stays **exact** for any
  surrogate quality (the metric is a preconditioner, not part of the target;
  the Metropolis correction with the exact density is the corrector — no delayed
  acceptance is needed). Importance weighting (`weights`, `1/Q`) and the
  closed-form-or-simulated Fisher slot (`fisher`) are exposed for the general
  case; `gdpar_geom_reservoir()` collects phase-one sites from a warmup run.
* **Implicit leapfrog hardened.** A numerical failure of either fixed-point
  solve (a step too large for the local curvature: a non-finite iterate, a
  Hessian that overflows, a non-positive-definite metric) is now caught and
  reported as non-convergence so the proposal is rejected, rather than crashing
  the run. The Euclidean and analytic-Fisher branches are unchanged
  (bit-identical, goldens intact).
* **Validated in three rigour layers.** Algebraic: the log-Cholesky
  parametrisation round-trips and both its forward and inverse differentials are
  exact; the closed-form metric derivative matches finite differences; the
  surrogate interpolates the supplied Fisher at the reservoir sites. Numerical:
  the implicit leapfrog with the learned metric is reversible to machine
  precision and conserves energy without drift. Methodological: the metric
  degrades continuously to SoftAbs far from the reservoir and the predictive
  standard deviation flags out-of-distribution positions. Statistical (gated):
  on Neal's funnel the learned metric reproduces the Capa 1 funnel exploration
  (E-BFMI lifted well above the Euclidean metric, recovering `Var(v)`).

## Bloque RG — RG.3 Riemannian metric + implicit generalised leapfrog (Session B9.24, 2026-06-02)

* **RG.3: position-dependent (Riemannian) geometry, the funnel remedy
  (`gdpar_geom_metric_riemannian()`).** The first level-3 metric of the Block RG
  hierarchy, extending the pluggable interface of the RG.2 engine with
  `position_dependent = TRUE`. Two curvature sources: the **expected Fisher
  information** (`curvature = "fisher"`, the natural metric of the statistical
  manifold, positive-definite by construction and supplied as a function), the
  primary maximally robust choice; and the **SoftAbs** regularisation of the
  observed Hessian (`curvature = "softabs"`, Betancourt 2013, eigenvalues mapped
  to `lambda * coth(alpha * lambda)`), fully general and used as the cold-start
  and extrapolation fallback. The metric exposes the spatial derivatives
  `dmass(theta)` (the SoftAbs derivative via the Daleckii–Krein formula; the
  Fisher derivative analytic or finite-differenced).
* **Generalised implicit leapfrog (Girolami & Calderhead 2011).** The
  position-dependent integrator the Euclidean engine previously aborted on: two
  fixed-point sub-steps (in the momentum, then the position) bracket an explicit
  momentum sub-step, exactly time-reversible and volume-preserving up to the
  fixed-point tolerance; non-converged solves are counted divergent and
  rejected. `gdpar_geom_hmc()` now records the energy trace and the **E-BFMI**
  (energy Bayesian fraction of missing information) and routes position-dependent
  metrics through this integrator automatically.
* **Validated in three rigour layers.** Deterministic: the implicit leapfrog is
  reversible to machine precision and conserves energy without drift; algebraic:
  the metric derivative matches finite differences and the kinetic gradient
  matches finite differences of the kinetic energy; statistical (gated): on
  Neal's funnel the Riemannian sampler lifts E-BFMI from about 0.08 to about 1.0,
  reaches deep into the neck and recovers `Var(v)`, where the Euclidean metric
  stalls; cmdstan-gated: Stan's `$hessian` matches the analytic Hessian of the
  suite targets and yields the same SoftAbs metric. The Riemannian sampler stays
  exact for any metric (the metric is a preconditioner, not part of the target;
  the Metropolis correction with the exact log-density is the corrector — no
  delayed acceptance is needed for the metric). Opt-in: the default fit path is
  untouched, goldens stay bit-identical.

## Bloque RG — RG.1.c calibration + RG.2 R-native geometric engine (Session B9.23, 2026-06-02)

* **RG.2: R-native geometric sampling engine (`gdpar_geom_hmc()`,
  `gdpar_geom_target()`, `gdpar_geom_metric_euclidean()`).** A pure-R
  Hamiltonian integrator that delegates the log-density, gradient and Hessian to
  a compiled backend — a cmdstan model built with
  `compile_model_methods = TRUE` (exposing `$log_prob` / `$grad_log_prob` /
  `$hessian`) or an R closure such as the dual targets of
  `gdpar_geometry_suite()`. This is the canonical decision-A scaffolding (D88)
  for the Block RG geometry hierarchy: pluggable metric / kinetic-energy /
  symplectic-integrator interfaces with the **Euclidean level** delivered now
  (constant metric, standard leapfrog, static HMC with a Metropolis
  correction). The Riemannian (Fisher / SoftAbs, RG.3), Finsler / relativistic
  and sub-Riemannian levels extend the same interfaces. Validated in three
  rigour layers: a **deterministic** check against the exact Hamiltonian flow of
  a Gaussian (energy conservation, reversibility, second-order accuracy), a
  **statistical** check that HMC recovers the easy targets' moments, and a
  cmdstan-gated **algebraic** check that the R closure gradient matches Stan's
  `$grad_log_prob` and that a cmdstan-backed target drives the same trajectory.
  Opt-in and standalone: the default fit path is untouched, goldens stay
  bit-identical.
* **RG.1.c: full calibration of the posterior-geometry diagnostic.** The
  "8/8 diagonal of one replica" smoke of RG.1 is replaced by **measured error
  rates** over a difficulty × pilot-budget × replica grid, with a minimax-robust
  adaptive replica allocation (replicas spent only where the Wilson interval for
  the correct rate is unresolved) and an out-of-sample (held-out) threshold
  calibration. Design-uniform accuracy rose from 0.78 to 0.91 and held-out
  balanced accuracy from 0.60 to 0.89. Reported honestly in
  `inst/benchmarks/results/block_rg_calibration.md`: six of the eight pathology
  classes are recovered at 0.93–1.00, while **funnel and heavy tails remain
  mutually confusable** (recall ~0.6–0.7) and the deepest funnels and heaviest
  tails are the hardest cases — gdpar's diagnostic can and does miss them.
* **Re-calibrated `gdpar_geometry_thresholds()` defaults.** Eight thresholds
  moved to catch the signals as the *short* diagnostic pilots actually attenuate
  them (e.g. `condition_high` 50→12, `heavy_kurtosis_high` 3→1.8,
  `boundary_prox_high` 0.10→0.02, `nslope_grows` 0.30→0.80). The calibration is
  against the synthetic suite, not a claim of optimality on real posteriors; the
  helper stays exposed for re-calibration.

## Bloque RG — Sub-phase RG.1: posterior-geometry diagnostic + synthetic suite (Session B9.22, 2026-06-02)

* **New opt-in forensic capability: `gdpar_geometry_diagnostic()`.** Probes the
  geometry of a posterior with cheap Hamiltonian pilots and classifies the
  sampling pathology, localises the culprit parameter(s), estimates the
  difficulty-vs-n behaviour, recommends the geometry level that remedies it,
  and estimates cost. It is the diagnostic half of the Block RG capability
  (geometric robustness of sampling), motivated by the eBird Tweedie count
  coordinate that the no-U-turn sampler cannot traverse (session B9.21). The
  default fit path is untouched: this is a standalone, opt-in tool, so the
  goldens stay bit-identical.
* **Size-invariant signals only.** The classifier reads the divergence rate,
  the minimum energy fraction of missing information (E-BFMI), the tree-depth
  saturation rate, the posterior condition number, and the adapted-step-to-scale
  ratio. It deliberately does **not** use R-hat or the effective sample size as
  decision signals: on short pilots those were the false positives of sessions
  B9.20/B9.21. A **difficulty-vs-n curve** (pilots at growing sizes) separates
  quasi-determinism (difficulty grows with n) from structural pathology
  (constant), which is the 28-GB nuance: a big dataset is slow without being
  geometrically hard.
* **New `gdpar_geometry_suite()`: eight synthetic geometries of known
  difficulty.** A falsifiable calibration backbone covering the full pathology
  taxonomy (isotropic control, anisotropic, Neal's funnel, heavy tails,
  quasi-deterministic canyon, multimodal, boundary-pegged, flat direction).
  Each target is dual (a Stan program plus an R log-density closure with an
  analytic gradient on the unconstrained scale), so the same geometry can be
  exercised by the cmdstan pilots and, later, by the R-native geometric engine
  of RG.2. Analytic gradients are cross-checked against finite differences to
  better than 1e-9.
* **New `gdpar_geometry_thresholds()`** exposes the classifier thresholds as
  tunable data for the calibration sweep.
* **Smoke validation: 8/8 on the single-replicate confusion diagonal.** Run
  end-to-end through real cmdstan pilots, the diagnostic classified all eight
  synthetic geometries correctly and routed each to its taxonomy remedy
  (funnel to Riemannian, heavy tails to Finsler/relativistic, quasi-determinism
  to sub-Riemannian, multimodality to tempering, the flat direction to
  reparametrise/eliminate). This is the diagonal on a single seed; the full
  error-rate calibration over replicates and a threshold sweep is the deferred
  RG.1.c step.
* **No default-path or API regressions.** Purely additive surface (two exported
  functions, a thresholds helper, two `print` methods). No existing fit path,
  golden, or vignette changed.

## Bloque 9 — Sub-bloque 9.1 closed: Stan-bound adversarial re-validation (Session B9.9, 2026-05-29)

* **Internal synthetic adversarial re-validation complete (Criterion 1).**
  The Stan-bound layers deferred from B9.7/B9.8 were executed and pass with
  no drift and no failures (`FAIL 0 | WARN 8 | SKIP 0 | PASS 76`): the
  statistical fit-fuzz invariances (INV-C: permutation invariance and
  Lipschitz continuity of the EB location estimate), the 12-family
  outcome-outlier stress (STR-C: every canonical family returns finite
  estimates or aborts with a canonical error, never a silent NaN/Inf), and
  the full bit-exact compare-path.
* **No bit-exact drift across 22 compare-able goldens.** Re-fitting the FB
  Path C K×p goldens (beta/gamma/student_t, seeds 91002-91004), the 5
  scalar EB 8.6.B goldens, and the 8.3.9 Tier1 (12 configs) + Tier2
  (Tweedie K=3, B-spline p=2) goldens all reproduced their frozen draws
  bit-exactly. The cluster-K4 Stan unification preserved reproducibility;
  no reseed was needed.
* **Candidate refinement noted (D77, deferred).** The first substantive
  run of the gated stress layer surfaced that `compute_diagnostics()`
  computes `max(rhat, na.rm = TRUE)` over a vector that can be entirely
  `NA` under degenerate outlier-laden data, yielding `-Inf` plus a warning
  for the reported max R-hat. The point estimate stays finite (the contract
  holds); the all-`NA` guard is a documented candidate fix, confined to the
  gated stress path, and does not affect the default suite or `R CMD check`.
* **No source/API changes this session.** B9.9 is validation-only; the
  package surface is unchanged. See
  `inst/benchmarks/results/block9_internal.md` §11 for the full report.

## Bloque 9 — Sub-bloque 9.1 dependence-robust inference (Axis 2) + D74 fix (Session B9.8, 2026-05-28)

* **Non-finite outcome guard hardened (D74).** The six outcome-validation
  sites in `gdpar()` and `gdpar_eb()` now reject `Inf` / `-Inf` (as well as
  `NA` / `NaN`) outcomes pre-fit with a `gdpar_input_error`. Previously an
  infinite outcome slipped past the `is.na()` guard (`is.na(Inf)` is
  `FALSE`) and surfaced late as an opaque `gdpar_eb_numerical_error`. The
  broadened, type-aware predicate (`!is.finite()` for numeric outcomes,
  `is.na()` otherwise) preserves behaviour on non-numeric outcomes.
* **New `gdpar_dependence_diagnostic()` (D75).** Quantifies serial
  dependence in the residuals of a scalar Empirical-Bayes fit (lag-1
  autocorrelation, Durbin-Watson, Ljung-Box). It turns the iid-violation
  risk from an invisible theoretical hazard into a measured quantity and
  gates the remedy below. Scalar EB path only this session; the spatial
  analogue (Moran's I) and the full-Bayes path are deferred.
* **New `gdpar_dependence_robust()` (D75).** Dependence-robust standard
  errors and percentile intervals via a temporal moving / circular block
  bootstrap that refits the model on contiguous-block resamples of the
  data. This is the working-independence + robust-variance stance of
  Liang & Zeger (1986): the point estimates are unchanged (consistent
  under a correct mean structure, not efficient), only the reported
  uncertainty is made robust to temporal dependence. The default block
  length follows the n^(1/3) rate (Künsch 1989); the Politis-White (2004)
  data-driven constant is deferred.
* **Honest scope, stated plainly.** gdpar does NOT model residual
  correlated noise or spatial random effects (Axis 1); that is deferred to
  a future block, evidence-gated by the external validation 9.2. The
  package's posterior uncertainty (EB and FB alike) is not robust to data
  dependence under the conditional-independence likelihood; the
  pre-existing pre-flight "block bootstrap" is an MCMC-autocorrelation
  device for the CP/NCP info-ratio test, not a data-dependence-robust SE
  estimator.
* **New test file `test-block9_dependence.R`.** Stan-free resampler algebra
  and API guards run by default (including a bit-exact resampler golden,
  D76); end-to-end fit-based behaviour (AR(1) detection, robust-SE table,
  seed-determinism) is gated by `GDPAR_RUN_BLOCK9_DEP_FITS`. The
  HMC-dependent bootstrap SEs are not frozen bit-exact (HMC is not
  reproducible across Stan toolchains).
* **`test-block9_stress.R` Section A** gains STR-A4 (an infinite outcome
  aborts pre-fit with `gdpar_input_error`, the D74 contract).

## Bloque 9 — Sub-bloque 9.1 internal adversarial re-validation, layers 1-3 + drift smoke (Session B9.7, 2026-05-28)

* **Decision C resolved to C.iv (layered defence with drift-gated
  reseed; canonized as D73).** The Path C fuzz protocol for the internal
  synthetic adversarial re-validation (Charter §2.2) combines a
  structural adversarial layer, a property-based invariance layer, and a
  bit-exact regression layer used as a drift detector; the full reseed
  (C.ii) is conditional on detected drift. Consulted batched at session
  start under upfront-consultations with a six-dimensional analysis.
* **New test file `test-block9_adversarial.R` (Layer 1; 11 test_that).**
  Structural adversarial fuzz over the canonical architecture: aliased
  grouping / collinear designs (C4, C4-bis), duplicate
  `prior_canonical_kind` D-ID violations, structurally degenerate priors,
  and latent HOM stratification. All aborts are pre-fit
  (`gdpar_input_error` / `gdpar_identifiability_error`) or `passed =
  FALSE` verdicts; no Stan compilation. Documents two validation
  boundaries: prior *content* (negative scale) is a Stan-side reject, and
  HOM is a fit-time (not pre-fit) condition.
* **New test file `test-block9_invariances.R` (Layer 2; 6 test_that, 4
  default + 2 gated).** Deterministic equivariances of the
  identifiability eigenstructure: row permutation and covariate
  rescaling / sign flip leave the verdict and conditioning invariant
  (exact up to floating-point summation order). Statistical fit-based
  invariances (permutation invariance, noise-injection continuity) are
  gated by `GDPAR_RUN_BLOCK9_FIT_FUZZ` (nocturnal B9.8 sub-unit).
* **New test file `test-block9_stress.R` (Layer 3; 6 test_that, 5
  default + 1 gated).** NA / missing robustness (NA in AMM covariate, NA
  and NaN in outcome all abort with `gdpar_input_error`) and
  extreme-but-finite covariate outliers (99.9% / 0.1% quantiles, 1e6
  scale) keep the conditioning finite. The 12-family outcome-outlier
  stress sweep is gated by `GDPAR_RUN_BLOCK9_STRESS_FITS` (B9.8).
* **Bit-exact drift smoke (Layer 3 / 4).** The freshest cluster-K4
  golden (`golden_fb_KxP_gaussian_polynomial_K2_p2`, seed 91001) was
  reproduced in isolation and matched its frozen RDS bit-exactly
  (400 × 815 draws). No drift: the cluster-K4 Stan unification preserved
  bit-exactness by construction, so the conditional full reseed (C.ii) is
  not triggered.
* **Canonical report `inst/benchmarks/results/block9_internal.md`**
  opened (criterion 1 of the Block 9 closure; updated across B9.7-B9.9).
* New default-running tests: 121 PASS / 0 FAIL; 3 gated SKIP. No package
  R/Stan source was modified in B9.7 (tests + report + docs only); two
  candidate robustness refinements (non-finite outcome pre-validation;
  formula-RHS NA handling) are documented as D74 deferred to consultation.

## Bloque 9 — Sub-bloque 9.3.d D69 substantive closure: FB KxP fit harness + Tier 1 goldens (Session B9.6, 2026-05-27)

* **Deuda colateral D69 cerrada operativamente** (canonizada en B9.5
  bajo decision N.a, ejecutada en B9.6): nuevo driver interno
  `.gdpar_fb_KxP_fit()` (`R/stan_codegen.R`) que compone
  `.build_amm_design_KxP()` + el assembler unificado
  `.assemble_stan_data_KxP(path = "FB")` (M.iv B9.6) + cmdstanr
  sampling sobre la piece `amm_canonical_pmulti_KxP.stan` (J.iv.A
  B9.5). Bypassa el guard publico `.gdpar_guard_multiparam_multivariate`
  (que sigue activo para `gdpar()` user-facing); init conservador
  `init = 0.5` para evitar rejecting initial values con las familias
  Path B (beta / gamma / student_t cerca de los bordes del link).
* **Assembler unificado M.iv** (`R/eb.R`): la signature de
  `.assemble_stan_data_KxP()` ahora acepta `path = c("EB", "FB")` y
  `cp_W = FALSE`. El branch EB preserva byte-identico el comportamiento
  de Sub-fase 8.6.D (use_W = 0 hardcoded por D39 + stan_id en {1, 3}
  por D40'); el branch FB lifta ambas restricciones y emite los
  campos W esperados por la KxP piece (`X`, `dim_W`, `d`,
  `W_per_kj_dim`, `W_type_id`, `W_n_knots_full`, `W_knots_full`,
  `W_degree`). Backward-compatible por default `path = "EB"`.
* **Bootstrap script real**: `data-raw/bootstrap_fb_goldens_KxP.R`
  reemplaza los 4 fit closures stub (D69 canonizada en B9.5) por
  closures reales que invocan `.gdpar_fb_KxP_fit()` con `iter_warmup =
  200`, `iter_sampling = 200`, `chains = 2`, sobre datos sinteticos
  pequenos (n = 80-100) por seed canonical:
  91001 gaussian K=2 p=2 / 91002 beta K=2 p=2 / 91003 gamma K=2 p=3 /
  91004 student_t K=3 p=2. Persiste un snapshot por golden con
  `draws_matrix` + `stan_data` + canonical metadata bajo
  `tests/testthat/data/golden_fb_KxP_<family>_polynomial_K<K>_p<p>.rds`
  y appendea 4 filas al manifest CSV con `sub_phase = "9.3.d"` +
  `status = "fitable"`.
* **Tests B9.6**: 7 nuevos `test_that` blocks en
  `tests/testthat/test-stan_codegen_canonical.R` (frozen md5 de los
  6 EB pieces relocadas + CANONICAL_HELPERS placeholder sanity +
  legacy top-level removidos + KxP EB intactos + dispatcher render
  sin leftover placeholders + assembler EB-side preservado +
  assembler FB-side acepta Path B family set + driver guards).
  Section 4 nueva en `tests/testthat/test-fb_goldens_9_3_d_KxP.R`:
  fitable RDS structure + draws_matrix has theta_ref_kp + stan_data
  carries FB extension fields (M.iv).

## Bloque 9 — Sub-bloque 9.3.a colateral L.iv.A.2: 6 EB-side templates cascade to canonical pieces (Session B9.6, 2026-05-27)

* **Cascade L heredada cerrada operativamente** (canonizada en B9.5
  bajo decision L.iv.A para `amm_distrib_K`, ejecutada en B9.6 bajo
  O.a para los 6 EB-side templates restantes con helpers inline). Los
  6 EB-side templates con inline `bspline_basis_eval` +
  `apply_W_basis_diff` fueron relocados a
  `inst/stan/_canonical_pieces/` con el placeholder
  `// {{CANONICAL_HELPERS}}` reutilizando la infraestructura G.iv de
  B9.4:
  * `amm_canonical_eb_marginal.stan` (md5 7a279bf087a0fd0a95d35a490814673f, 6442 bytes)
  * `amm_canonical_eb_marginal_multi.stan` (md5 9250144e95df5fd409cb5cc9df3f3730, 13661 bytes)
  * `amm_canonical_eb_marginal_K.stan` (md5 dca0b8282e4c617d259e8f17d1f95be5, 36709 bytes)
  * `amm_canonical_eb_conditional.stan` (md5 8321f9d7d99d0c702045a2886cd1f98c, 6984 bytes)
  * `amm_canonical_eb_conditional_multi.stan` (md5 186b6ac9f9810584ea21522974fe641d, 13878 bytes)
  * `amm_canonical_eb_conditional_K.stan` (md5 169274eca2a64e2a0d5a060e2d0c3033, 36250 bytes)
* **Bit-exact Stan semantics preservadas** token-wise (6/6); solo
  difieren comments del helpers block, unified al inyectar la piece
  canonica `amm_canonical_helpers.stan` via el dispatcher.
* **Callers actualizados**: `generate_stan_code_multi()` +
  `generate_stan_code_K()` extendidos con dos `switch` entries cada
  uno traduciendo el legacy nombre al canonical (mirror del patron
  L.iv.A B9.5). `.gdpar_eb_render_template()` extendido con
  translation + helpers substitution (mirror del dispatcher canonical).
* **6 legacy top-level eliminados** ([[feedback_depuracion_obsoletos]]):
  `amm_eb_marginal{,_multi,_K}.stan` + `amm_eb_conditional{,_multi,_K}.stan`.
  Los 2 KxP EB templates (`amm_eb_marginal_KxP.stan` +
  `amm_eb_conditional_KxP.stan`) permanecen en `inst/stan/` root
  intactos: NO tienen helpers inline porque la canonization 8.6.D
  hardcodea `use_W = 0` per D39.
* **Cluster K4 cerrado completamente** (FB + EB cascade + fit
  harness): 4/4 deudas sustantivamente cerradas + L cascade
  propagada + D69 operativamente cerrada.



* **New canonical piece**
  `inst/stan/_canonical_pieces/amm_canonical_pmulti_KxP.stan`
  (~740 lines): Full-Bayes multi-parametric multivariate template for
  the K >= 2 K-individual slots crossed with p >= 2 AMM coordinates
  regime, the bit-exact FB counterpart of the EB Path C templates
  canonized in Sub-fase 8.6.D. Mirrors `amm_eb_marginal_KxP.stan`
  extended with the globally shared W modulating component and the
  full Path B family set
  `{1 gaussian, 3 neg_binomial_2, 5 beta, 6 gamma,
    7 lognormal_loc_scale, 8 student_t, 9 tweedie, 10 zip, 11 zinb,
    12 hurdle_poisson, 13 hurdle_neg_binomial_2}`. Canonized under
  decision I.iv lateral (DESIGN_9_3_D_PATH_C.md APERTURA cerrada
  B9.4) + J.iv.A piece arquitectonica (sub-decision 3.2 dedicated
  piece with canonical helpers reuse via `// {{CANONICAL_HELPERS}}`
  placeholder consumed by the dispatcher per the G.iv pattern of
  B9.4).
* **Dispatcher extension**: `.gdpar_emit_canonical_stan()`
  (`R/stan_codegen.R`) now accepts
  `p_class = "pmulti_KxP"` (plus `"distrib_K"` from the L.iv.A
  cascade below). The KxP branch consumes the canonical priors +
  `{{THETA_REF_PRIOR_BLOCK}}` (per-family anchor block via the new
  helper `.gdpar_build_theta_ref_prior_block_KxP()`; Tweedie K=3
  carries the bounded p slot per coord via the iterated slice form)
  + `{{W_SCALE}}` / `{{W_PRIOR}}` for the CP / NCP toggle. The KxP
  piece emits NCP-only a/b blocks inline (no
  `DATA_CP_A_PER_K_DECL` / `TP_A_BLOCK` / `MODEL_A_BLOCK` overrides
  for B9.5 atomic scope; per-slot per-coord CP / mixed variants are
  deuda for B9.6+).
* **Thin wrapper**: new internal `generate_stan_code_KxP_FB()`
  preserves the canonical source-of-truth pattern of B9.3; the
  user-facing guard `.gdpar_guard_multiparam_multivariate` still
  aborts `gdpar()` on K > 1 + p > 1, deferring the public surface
  lift to a future sub-phase.
* **Reserved canonical seeds** for the bootstrap goldens roster
  (DESIGN_9_3_D_PATH_C.md §4): `91001..91004` cover the K.c roster
  minimum (gaussian K=2 p=2, beta K=2 p=2, gamma K=2 p=3,
  student_t K=3 p=2). `91005..91099` stay reserved for the
  expanded roster diferido a B9.6.
* **Bootstrap script** `data-raw/bootstrap_fb_goldens_KxP.R`
  (env-var gated on `GDPAR_BOOTSTRAP_FB_GOLDENS` +
  `GDPAR_BOOTSTRAP_FB_GOLDENS_KxP`, mirroring the 8.6.D pattern):
  records 4 metadata-only manifest rows for the reserved seeds. The
  actual cmdstanr sampling is deuda D69 canonized in B9.5 +
  deferred to B9.6 under `[[feedback_token_economy_session_planning]]`;
  the B9.5 cycle closed the codegen + dispatcher + piece architecture
  with bit-exact tests over the canonical Stan source (frozen md5
  invariants).
* **New roster test**
  `tests/testthat/test-fb_goldens_9_3_d_KxP.R`: 8 codegen tests +
  3 RDS-presence tests (skip_if RDS not present + skip_if manifest
  not present, the canonical 8.6.D convention). Always-on codegen
  tests verify the family-specific dispatch + the seed range
  reservation + the bootstrap script env-gating + the determinism
  of the wrapper.
* **Test extensions**:
  `tests/testthat/test-stan_codegen_canonical.R` adds frozen md5
  invariants for the new KxP piece + the relocated distrib_K piece +
  the dispatcher OUTPUT for gaussian / tweedie + helpers presence
  sanity + cmdstanr parse-check (skip_if cmdstan unavailable).
* **Caveat heredado Path B logit-strict bajo Path C** registered in
  the new piece header and in DESIGN_9_3_D_PATH_C.md §5: the K.c
  B9.5 roster (gaussian + beta + gamma + student_t) avoids
  logit-strict puro; roster expansion to {10 ZIP, 12 hurdle_poisson}
  (link logit in pi) reserved for B9.6 with skip_if granular if
  emergent numerical instability appears.

## Bloque 9 — Sub-bloque 9.3.a colateral L.iv.A: amm_distrib_K cascade to canonical pieces (Session B9.5, 2026-05-27)

* **Relocated canonical piece**
  `inst/stan/_canonical_pieces/amm_canonical_distrib_K.stan`
  (851 lines): byte-relocated copy of the legacy
  `inst/stan/amm_distrib_K.stan` with the inline
  `bspline_basis_eval` + `apply_W_basis_diff` definitions replaced
  by the `// {{CANONICAL_HELPERS}}` placeholder consumed by the
  dispatcher per the G.iv pattern of B9.4. The Tweedie family
  helpers (`tweedie_log_W_series`, `tweedie_log_f_series`,
  `tweedie_log_f_saddlepoint`, `tweedie_lpdf`, `tweedie_rng`) and
  `apply_inv_link_by_id` stay inline because they are K-specific
  (not shared across templates).
* **Dispatcher extension**: `.gdpar_emit_canonical_stan()` accepts
  `p_class = "distrib_K"` with the full K-style placeholder set
  (`THETA_REF_PRIOR_BLOCK` + `DATA_CP_A_PER_K_DECL` + `TP_A_BLOCK` +
  `MODEL_A_BLOCK` + canonical priors + W toggles). The
  `template_name` override supports the EB-side K templates
  (`amm_eb_marginal_K.stan`, `amm_eb_conditional_K.stan` per D34 of
  Sub-fase 8.6.C); those EB templates retain inline helper copies
  pending the EB cascade in B9.6 (deuda L heredada bajo L.iv.A
  symmetry).
* **`generate_stan_code_K()` refactor**: now a thin wrapper that
  delegates to the canonical dispatcher with `p_class = "distrib_K"`,
  preserving signatures and tests bit-for-bit. The legacy
  `template_name = "amm_distrib_K.stan"` default is translated
  internally to `"amm_canonical_distrib_K.stan"`; the EB-side K
  templates pass through unchanged.
* **Legacy deletion**: `inst/stan/amm_distrib_K.stan` removed per
  `[[feedback_depuracion_obsoletos]]`.
* **Bit-exact preservation**: the Stan semantics of
  `generate_stan_code_K(prior, family = ...)` are preserved bit-exact
  across L.iv.A (only the unified helper-block comments differ vs the
  legacy template). The 14 K=2 / K=3 goldens in
  `tests/testthat/data/golden_K*.rds` remain valid by construction
  (they fit data, not source hashes); cmdstanr cache recompiles
  once on first call post-B9.5 and is then warm.
* **Caveat L cascade**: the 6 EB-side templates with inline
  helpers (`amm_eb_marginal*.stan`, `amm_eb_conditional*.stan` for
  single / multi / K variants) retain inline copies pending B9.6
  cascade (canonized as deuda colateral L heredada in
  CHARTER_BLOQUE_9.md §2.1).

## Bloque 9 — Sub-bloque 9.3.c: gdpar_ksd_joint() helper for joint KSD between EB and FB posteriors (Session B9.4, 2026-05-27)

* **New exported helper `gdpar_ksd_joint(eb_fit, fb_fit, ...)`**
  (`R/ksd_joint.R`): operationalizes the open question documented in
  the Roxygen of `gdpar_compare_eb_fb()` that the marginal
  total-variation distance is only a coarse proxy and the joint
  posterior discrepancy deserves a density-free spectral metric.
  Returns an object of class `gdpar_ksd_joint` with the KSD
  V-statistic, kernel/bandwidth configuration, target empirical
  Gaussian (mean + covariance), and call. Canonized under decision
  H.iv lateral: IMQ base kernel default (Gorham-Mackey 2017,
  dim-independent rate for log-concave targets) + median heuristic
  bandwidth + ESS-weighted variant via `posterior::ess_basic`. RBF
  base kernel is provided as a textbook alternative.
* **Target choice (B9.4 iteration)**: empirical Gaussian Laplace
  approximation of the FB posterior over the common xi-variables;
  closed-form Stein kernel via Sigma_hat^-1. The full-KSD variant
  against the actual FB target via cmdstanr's `grad_log_prob()` is
  documented as a Block 9.x extension.
* **S3 methods**: `print.gdpar_ksd_joint`, `summary.gdpar_ksd_joint`,
  `print.summary.gdpar_ksd_joint`.
* **Hot fix during B9.4 suite verification**: `.gdpar_ksd_safe_min_ess`
  was updated to use the canonical `posterior::summarise_draws()` API
  (was calling `posterior::ess_basic()` directly on a `draws_matrix`,
  which is method-dispatched ambiguously across posterior versions and
  could return either a numeric scalar or a tibble). Resolves the two
  ESS-weighted-related errors observed in the full-suite run.
* **Vignette section**: new Section 11 in
  `v07b_eb_multivariate.Rmd` ("Joint Kernel Stein Discrepancy: the
  `gdpar_ksd_joint()` Helper") with subsections on the Stein
  operator, the empirical Gaussian target, the base kernel and
  bandwidth, the ESS-weighted variant, and the complementarity
  between `gdpar_compare_eb_fb` (marginal TV) and `gdpar_ksd_joint`
  (joint KSD). Existing Section 11 (Summary) renumbered to Section
  12.
* **Tests** (`tests/testthat/test-ksd_joint.R`): 14 test_that blocks
  / 45 assertions covering input validation (5), happy path and
  class invariants (3), statistical layer (KSD small under matched
  Gaussians + KSD detects shifted Gaussians) (2), ESS-weighted
  variant (2), S3 methods (2). Mock fits avoid Stan dependence.

## Bloque 9 — Sub-bloque 9.3.a colateral: deuda f closed via dedicated canonical helpers piece (Session B9.4, 2026-05-27)

* **Canonical helpers piece `inst/stan/_canonical_pieces/amm_canonical_helpers.stan`**
  (new file): dedicated Stan source FRAGMENT containing the two
  canonical helpers `bspline_basis_eval` and `apply_W_basis_diff`.
  Closes deuda f of 9.3.a (helpers duplication across the K=1 FB
  pieces) under canonized decision G.iv lateral (R-side textual
  substitution). The piece has NO surrounding `functions { }` block;
  it is intentionally not standalone-parseable by `cmdstanr::stanc`
  and is inserted by `.gdpar_emit_canonical_stan()` at the
  `// {{CANONICAL_HELPERS}}` placeholder inside the body piece's
  `functions { }` block.
* **Body pieces shrunk**: `amm_canonical_p1.stan` (13878 → 10722
  bytes, md5 `f899811e...` → `730c93f5...`) and
  `amm_canonical_pmulti.stan` (16028 → 13764 bytes, md5
  `112a4e68...` → `01fc4b04...`) now contain only the
  `// {{CANONICAL_HELPERS}}` placeholder inside their `functions { }`
  block; helper definitions removed inline. Total bytes saved across
  the K=1 FB pieces: 2240 (~7.5%); future EB-side cascade in B9.5+
  will reuse the same helpers piece with additional savings.
* **Dispatcher extended**: `.gdpar_emit_canonical_stan()` reads the
  helpers piece and substitutes the placeholder via `gsub()` before
  applying the other prior / parametrization placeholders. The gate
  is `grepl("// {{CANONICAL_HELPERS}}", src, fixed = TRUE)`, so
  templates that do not contain the placeholder are unaffected
  (forward-compatible with the EB cascade and with any future piece
  that ships its own helpers).
* **Bit-exact OUTPUT preserved for the p1 path** (md5
  `8930c724...` for NCP and `1e6d7d28...` for MLE):
  the helpers in the dedicated piece are verbatim from the
  pre-G.iv `amm_canonical_p1.stan`. **Output md5 changes for the
  pmulti path** (`cf182996...` → `2ede60f6...`, +892 bytes) due to
  unified helper comments (pre-G.iv `amm_canonical_pmulti.stan`
  carried a brief 4-line block referencing the then-open deuda f;
  it is replaced by the detailed 16-line algorithmic comment block).
  **Stan semantics preserved bit-exact** (only comments differ);
  goldens unaffected (they fit data, not source hashes);
  `cmdstanr::cmdstan_model()` cache hash recompiles once on first
  call post-B9.4.
* **Tests** (`tests/testthat/test-stan_codegen_canonical.R`): 5 new
  test_that blocks (frozen md5 of helpers piece + frozen output
  md5 invariants for p1 NCP / p1 MLE / pmulti default + sanity
  check that helpers are inserted + defensive smoke). The 2 existing
  frozen-piece-md5 tests updated with the new post-G.iv hashes.

## Bloque 9 — Sub-bloque 9.3.a: unification of K=1 FB Stan templates via canonical R-side codegen dispatcher (Session B9.3, 2026-05-27)

* **Canonical dispatcher `.gdpar_emit_canonical_stan()`** (new internal
  helper, `R/stan_codegen.R`): single source-of-truth dispatcher for
  the (p, K=1) FB Stan templates. Replaces the duplicated
  template-read plus substitute logic that previously lived inline in
  `generate_stan_code()` (p = 1) and `generate_stan_code_multi()`
  (p >= 1). Reads a canonical piece selected by `spec$p_class` and
  applies the appropriate placeholder substitution. The dispatcher
  emits the same substituted Stan strings as the legacy paths for the
  same inputs.
* **Pieces relocated to `inst/stan/_canonical_pieces/`** (new
  directory): the legacy `inst/stan/amm_main.stan` and
  `inst/stan/amm_distrib_multi.stan` are renamed to
  `amm_canonical_p1.stan` and `amm_canonical_pmulti.stan` respectively
  and moved into `inst/stan/_canonical_pieces/`. The relocation is
  byte-identical (md5sum preserved); bit-exactness of the substituted
  Stan source is preserved by construction; no golden re-bootstrap
  required.
* **Legacy top-level templates removed**: `inst/stan/amm_main.stan`
  and `inst/stan/amm_distrib_multi.stan` are deleted from the package
  root. Backward-compat for callers of `generate_stan_code_multi()`
  that pass `template_name = "amm_distrib_multi.stan"` is preserved
  via wrapper-side translation to the canonical piece name; other
  `template_name` values (the EB-side multi templates pending the
  cascade in B9.5+: `amm_eb_marginal_multi.stan`,
  `amm_eb_conditional_multi.stan`) are served from
  `inst/stan/` root unchanged.
* **Public-internal wrappers thinned**: `generate_stan_code()` and
  `generate_stan_code_multi()` become thin wrappers (~10 lines each)
  that construct a canonical codegen spec and delegate to
  `.gdpar_emit_canonical_stan()`. Signatures preserved bit-for-bit;
  no caller migration required across the package. The K-individual
  generator `generate_stan_code_K()` (K >= 2, p = 1, multi-parametric
  via `amm_distrib_K.stan`) is out of scope for 9.3.a and is left
  untouched.
* **New tests** (`tests/testthat/test-stan_codegen_canonical.R`,
  9 tests): frozen md5 hashes of the two canonical pieces;
  determinism of the dispatcher across calls (both p1 and pmulti
  paths); equivalence of the public wrappers with the direct
  dispatcher call; input validation (invalid `p_class`, malformed
  `cp_a_per_k`); cmdstanr parse-check of the substituted output for
  both pieces (skipped when cmdstan is unavailable).
* **Updated tests** (`tests/testthat/test-stan_codegen_multi.R`): the
  legacy existence check for `amm_distrib_multi.stan` is canonized to
  verify the relocation (canonical pieces exist; legacy top-level
  templates are removed).
* **Decision B canonized** (B.iv lateral): codegen R-side as the
  canonical source of truth wins 6/6 on the multidimensional
  robustness analysis (conceptual + theoretical + statistical +
  computational with design-robustness + mathematical + methodological)
  against the alternative strategies of parallel introduction with
  guard (B.i), direct replacement (B.ii), and indefinite hybrid
  preservation (B.iii). Sub-decisions A1.iv (refactor in-place +
  canonical dispatcher with thin wrappers), A2 (pieces in
  `_canonical_pieces/`), A3.1 (snapshot frozen + md5 hash dual
  verification), B4 (FB cluster first, EB cascade in B9.5+), B5
  (deletion during B9.3 once tests pass), B6 (9.3.c intercalable in
  B9.4), C7.2 (deuda f helpers dedup deferred to B9.4 under the
  token-economy principle), D8 (rename without lifecycle since
  functions are internal), E9 (Roxygen + this NEWS entry; vignette
  section in `v07_path_1.Rmd`), F10 (closure criterion under rigor
  de tres capas), F11 (atomic sub-unit pacing per
  `feedback_token_economy_session_planning`) canonized concurrently.
* **Sub-deuda f (helpers `bspline_basis_eval` and `apply_W_basis_diff`
  duplicated in pieces)** remains open as planned for B9.4 under the
  token-economy session-planning principle; `#include` mechanism to
  factor the helpers into `_canonical_pieces/amm_canonical_helpers.stanfunctions`
  will be canonized in the next session.

## Bloque 9 — Sub-bloque 9.3.b: O5\*-EBFB Hessian anti-fragility (Session B9.2, 2026-05-27)

* **Adaptive Levenberg-Marquardt ridge** (canon Section 2.8 of
  `CHARTER_SUBFASE_8_6.md`, component 2): the single-step ridge of
  Sub-phase 8.6.B is extended into an adaptive geometric loop bounded
  by `laplace_control$ridge_max_iter` (default 10) and multiplied by
  `laplace_control$ridge_grow_factor` (default 10.0) at each iteration
  until the post-ridge condition number is at or below
  `laplace_control$kappa_threshold` or the loop is exhausted.
* **`|det(cov)| < epsilon_lm` trigger** (new condition): the helper
  also triggers when the determinant of the empirical Laplace
  covariance is strictly below `laplace_control$epsilon_lm`
  (default `sqrt(.Machine$double.eps)`, approximately 1.5e-8) even
  when the matrix is technically positive-definite, matching the
  canonical Section 2.8 wording \emph{|det(H)| < epsilon_LM} of the
  ridge trigger.
* **`.gdpar_eb_lm_perturb()`** (new internal helper, `R/eb.R`):
  factored out from the per-Path inline ridge of Sub-phase 8.6.B/C/D,
  shared by both `.gdpar_eb_maximize_marginal()` (Paths A/B) and
  `.gdpar_eb_maximize_marginal_KxP()` (Path C, per-slot). Returns a
  list with `cov_perturbed`, `lambda_used`, `n_iter`, `kappa_post`
  and `status` (one of `"not_needed"`, `"converged"`, `"exhausted"`).
* **Path C now applies the ridge per-slot** (previously absent): each
  of the K slot covariances in the Path C K x p regime is perturbed
  independently with its own lambda, matching the per-slot canonical
  aggregation already established by D43 of Sub-phase 8.6.D.
* **New `diagnostics_numerical` fields**: `lm_n_iter`, `lm_status`,
  `kappa_post_ridge` (Paths A/B); `lm_lambda_per_slot`,
  `lm_n_iter_per_slot`, `lm_status_per_slot` (Path C). The pre-9.3.b
  `lm_perturbation` and `kappa_per_slot` fields are preserved.
* **Updated `gdpar_eb_numerical_error` messages** include the L-M
  status, lambda used and number of iterations for both Paths A/B and
  Path C, enabling downstream code to disambiguate `"exhausted"` from
  `"converged but still above threshold"`.
* **Tests**: 8 new unit-level tests in
  `tests/testthat/test-eb_numerical_robustness.R` exercising the new
  helper across well-conditioned, rank-deficient, small-determinant,
  exhausted-by-budget and 1x1 covariances; plus 4 new resolver tests
  for the three new `laplace_control` defaults and their validators.

## Sub-phase 8.6.E — EB-vs-FB comparator + operational vignette + aggregated 8.6 closure (2026-05-26)

* **`gdpar_compare_eb_fb(eb_fit, fb_fit, level = 0.95, tv_bins = 30L)`**
  (new function in `R/compare_eb_fb.R`): canonizes the operational
  EB-vs-FB comparison promised in v07 Section 11.1. Takes a
  `gdpar_eb_fit` and a `gdpar_fit` fitted on the same dataset and
  returns a `gdpar_eb_fb_comparison` object with three tables:
  * `theta_diff_table`: per-anchor cell comparison of
    $\widehat\theta_{\text{ref}}^{\text{EB}}$ vs the FB posterior mean
    $E_{\text{FB}}[\theta_{\text{ref}}]$, with the difference and the
    difference normalized by the FB standard error. Operationally
    verifies Theorem 7A / 7A\* (first-order asymptotic equivalence).
  * `tv_table`: marginal empirical total variation distance per common
    $\xi$ parameter via histogram plug-in. Theorem 7A predicts marginal
    TV $\to 0$ in probability as $n \to \infty$.
  * `coverage_table`: per anchor cell, EB credible-interval width
    (with inflation applied when `eb_correction = TRUE`), FB
    credible-interval width, and the ratio
    `width_eb / width_fb`. Operationally verifies the $O(n^{-1})$
    under-cover prediction of Proposition 7B / 7B\* matricial /
    7B\* tensorial.
* **Path coverage of the comparator**: handles all four EB regimes
  uniformly (base $K=1+p=1$, Path A $K=1+p>1$, Path B $K>1+p=1$,
  Path C $K>1+p>1$). Path C tables are keyed by (group, slot, coord);
  base / Path A / Path B by cell index. The K x p tensor inflation
  contributes per-cell inflation factors via the diagonal of each
  slot's correction block.
* **S3 methods**: `print.gdpar_eb_fb_comparison`,
  `summary.gdpar_eb_fb_comparison`,
  `print.summary.gdpar_eb_fb_comparison` in
  `R/compare_eb_fb_methods.R`. The summary aggregates the TV table
  and the coverage table into quartile / min / median / max / mean
  summaries with the full per-anchor diff table.
* **`vignettes/vop07_eb_workflow.Rmd`** (new operational vignette,
  ~400 lines): walks through `gdpar_eb()` end to end across the four
  path regimes, the numerical diagnostics ($\kappa(H)$, LM ridge,
  multi-start dispersion, `kappa_per_slot` and slot dispositions under
  Path C), the EB-vs-FB comparison via `gdpar_compare_eb_fb()`, and
  the troubleshooting recipes for ill-conditioned Hessians and
  conditional-HMC instability under logit-strict links (Path B
  caveat heredado de 8.6.C + Path C caveat heredado de 8.6.D).
  Chunks default to `eval = FALSE` because each fit compiles Stan
  models and takes minutes; re-enable per-chunk or globally to
  reproduce the runs locally.
* **`vignettes/v07_eb_vs_fb.Rmd`** Section 11.1 rewritten in
  **descriptive** voice apuntando a `gdpar_eb()` con las cuatro
  variantes de path Stan templates y a `gdpar_compare_eb_fb()`; se
  elimina el contrato textual prospectivo `eb_for_theta_ref = TRUE`
  (nunca implementado per Charter Decision 2.2). Section 13
  («Connections to Subsequent Blocks») extendida con dos viñetas
  apuntando a `v07b_eb_multivariate` (extensión teórica multivariada)
  y `vop07_eb_workflow` (operativa).
* **`tests/testthat/test-compare_eb_fb.R`** (new): unit tests for the
  comparator with 4 input-validation guards (unconditional) + 1
  end-to-end smoke gated by `GDPAR_RUN_STAN_SMOKE_EB=1` + 4
  internal-helper unit tests for `.gdpar_eb_fb_tv_table` /
  `.gdpar_eb_fb_extract_xi_draws`.
* **`NAMESPACE`**: exports `gdpar_compare_eb_fb` + 3 S3method
  registrations regenerated via `devtools::document()`.
* **`man/`**: `gdpar_compare_eb_fb.Rd` +
  `print.gdpar_eb_fb_comparison.Rd` +
  `summary.gdpar_eb_fb_comparison.Rd` +
  `print.summary.gdpar_eb_fb_comparison.Rd` regenerated.

This closure aggregates Sub-fase 8.6 entera (A theoretical extension
+ B base implementation + C multivariate partial paths + D full
Path C + E comparator and operational vignette). For the canonical
opening + canonization summary of Sub-fase 8.6 see
`CHARTER_SUBFASE_8_6.md`; for the closure summary across the five
sub-sub-fases see `HANDOFF_SUBFASE_8_6_CIERRE.md`.

## Sub-phase 8.6.D — Empirical Bayes full multivariate Path C ($K > 1 \wedge p > 1$) (2026-05-26)

* **Path C** ($K > 1 \wedge p > 1$) of `gdpar_eb()` is now operative.
  The path inherits two new Stan templates canonized in Sub-phase 8.6.D
  under **decision D36 = (alpha)** (Session 13a 2026-05-25):
  `inst/stan/amm_eb_marginal_KxP.stan` and
  `inst/stan/amm_eb_conditional_KxP.stan`. The canonical Path C
  architecture is the composition of Path A multivariate ($p > 1$,
  Theorem 7A* of v07b) with Path B multi-parametric ($K > 1$, Theorem
  7C* of v07b) under **decision D38'' = (h)** (Session 13b 2026-05-25):
  outcome `y[n, p]` matrix-column shared across the K slots, with
  per-slot per-coord linear predictors `eta_kp[k, i, j]` for the
  canonical Path B K-slot multi-parametric family.
* **Decision D37 = (i)** (tensor full $K \times p \times p$ correction):
  `.gdpar_eb_correction_tensor()` returns a 3D array
  $C^{\text{tensor}}_{g, \alpha}[k, , ] = \kappa(\alpha) \cdot
   \Sigma^{\text{marg}}_{\theta_{\text{ref}, k}}$
  retaining all cross-coordinate terms within each slot. The
  block-diagonal structure across slots is the canonical factorization
  assumption $\pi_\xi = \prod_k \pi_{\xi_k}$. Reduces algebraically to
  the Path A matrix correction at $K = 1$ and to the scalar Path B
  correction at $p = 1$. The `correction_tensor_constant` /
  `correction_tensor_dispositions` slots of `gdpar_eb_fit` report the
  tensor and per-slot dispositions (`"ok"`, `"non_finite"`,
  `"non_psd"`, `"missing"`, `"disabled"`).
* **Coverage initial under decision D40'** (Session 13b 2026-05-25):
  Path C supports `stan_id %in% c(1, 3)` per slot (Gaussian K=2,
  Negative Binomial K=2). The remaining Path B set
  $\{5, 6, 7, 8, 9, 10, 11, 12, 13\}$ (Beta / Gamma /
  Lognormal-loc-scale / Student-t / Tweedie / ZIP / ZINB /
  Hurdle-Poisson / Hurdle-NB) is deferred to a later iteration of 8.6.D
  with the explicit numerical caveat of Section 6.1 of the opening
  handoff (HMC condicional bajo plug-in EB cerca del borde de soporte
  logit/log links + warmup corto). The `.gdpar_eb_check_stan_id_for_path`
  guard rejects deferred families with a clear deferral message.
* **Decision D39** (Session 13b 2026-05-25): `use_W = 0` enforced in
  Path C first iteration; the modulating component W is registered as
  a follow-on debt against Block 9.x because its canonical
  factorization across the $K \times p$ cross structure has no clean
  canonization in v07b yet. The orchestrator upstream rejects input
  with `W != NULL` on any slot with `gdpar_unsupported_feature_error`.
* **Decision D41** (Session 13b 2026-05-25): `a` / `b` coefs are
  per-slot per-coord ragged in Path C. Each `(slot k, coord j)` pair
  has its own `J_a_per_kp[k, j]` and `J_b_per_kp[k, j]`; the design
  matrices `Z_a_kp[k, j]` and `Z_b_kp[k, j]` are padded matrices
  `matrix[n, J_a_max]` for maximum statistical generality without
  coord-wise restrictions within a slot.
* **Decision D43 = (a)** (Session 13c 2026-05-26): the Step (i)
  Laplace strategy under Path C is a single joint Laplace
  approximation over the full tensor of `theta_ref_kp`
  $[J_{\text{groups}} \times K \times p]$ anchors followed by
  block-diagonal extraction per slot post-hoc to feed
  `.gdpar_eb_correction_tensor()`. Minimum compile count (1 vs K
  separate Stan compiles for the per-slot alternative); preserves the
  EB-Path-B-Compound canonical estimator (Petrone-Rousseau-Scricciolo
  2014). The block-diagonal extraction averages over groups when
  $J > 1$ to match the per-slot canonical aggregation of v07b
  Section 5.1.
* **End-to-end Path C orchestrator**: `.gdpar_eb_run_KxP()` (R/eb.R,
  ~350 lines) drives the full pipeline: input validation +
  homogeneous-p enforcement, design building via
  `.build_amm_design_KxP()`, stan_data assembly via
  `.assemble_stan_data_KxP()`, conservative path-aware init via
  `.gdpar_eb_make_random_init_KxP()` (Path C-specific because the
  default cmdstanr unconstrained random init explodes the K x p joint
  via `exp(eta_kp[2])` overflow), Laplace via
  `.gdpar_eb_maximize_marginal_KxP()`, conditional HMC sampling, and
  tensor correction application. Returns a `gdpar_eb_fit` with
  `path = "eb_KxP"` and new slots `theta_ref_kp_hat`
  (3D array $[J, K, p]$), `theta_ref_kp_se`, `theta_ref_kp_cov_per_slot`
  (list of K matrices $p \times p$), `correction_tensor_constant`,
  `correction_tensor_dispositions`.
* **`print.gdpar_eb_fit` and `summary.gdpar_eb_fit`** extended with a
  Path C branch: report per-slot per-coord `theta_ref_kp_hat` and SE
  for group 1, per-slot kappa diagnostics, and per-slot dispositions.
  The summary tabulates the corrected credible intervals per `(group,
  slot, coord)` cell with the per-slot inflation factor derived from
  the diagonal of the correction tensor block.
* **`R/golden_helpers.R`**: `.gdpar_golden_roster_8_6_D()` canonizes 8
  configurations (4 fitable Gaussian / NB at $p \in \{2, 3\}$ + 4
  metadata-only guard Beta / Gamma / Student-t / ZIP at $p = 2$).
  `.gdpar_golden_eb_rds_name_8_6_D(family, regime)` produces the
  canonical filename
  `golden_eb_<family>_pathC_<regime>.rds` so the 8.6.B / 8.6.C / 8.6.D
  on-disk namespaces never collide.
* **`tests/testthat/test-eb_smoke_by_family.R`**: 8 new smokes appended
  (smokes 35-42). 4 Path C guard smokes (39-42) run unconditionally
  (verify `gdpar_unsupported_feature_error` on Beta / Gamma /
  Student-t / ZIP under $K, p > 1$). 4 Path C fitable smokes (35-38;
  Gaussian K=2 and NB K=2 crossed with $p \in \{2, 3\}$) are gated by
  `GDPAR_RUN_STAN_SMOKE_EB=1`.
* **`tests/testthat/test-eb_goldens_8_6_D.R`**: new structural test
  mirroring `test-eb_goldens_8_6_C.R`. Validates roster coverage,
  Path C invariants ($K > 1 \wedge p > 1$ in every entry), filename
  convention, correction tensor shape under valid PSD inputs and under
  injected non-PSD slots, `.gdpar_eb_check_stan_id_for_path` Path C
  acceptance ($\{1, 3\}$) and Path B preservation
  ($\{1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13\}$ at $p = 1$), guard class
  canonicity, and manifest CSV cross-validation.
* **`data-raw/bootstrap_eb_goldens_D.R`**: new bootstrap script gated
  by `GDPAR_BOOTSTRAP_EB_GOLDENS=1` + `GDPAR_BOOTSTRAP_EB_GOLDENS_D=1`.
  Bootstraps the 4 fitable Path C goldens with seed `20260526L` and 4
  metadata-only guard goldens; appends 8 rows to
  `tests/testthat/data/golden_manifest.csv` with `sub_phase = "8.6.D"`.

## Sub-phase 8.6.C — Empirical Bayes multivariate partial paths (2026-05-25)

* **Path A** ($K = 1$, $p > 1$) and **Path B** ($K > 1$, $p = 1$) of
  `gdpar_eb()` are now operative. The two paths inherit four new Stan
  templates canonized in Sub-phase 8.6.C under **decision D34**:
  `inst/stan/amm_eb_marginal_multi.stan` and
  `inst/stan/amm_eb_conditional_multi.stan` (Path A; bit-identical bodies
  to `amm_distrib_multi.stan` modulo the `theta_ref_data` move + EB
  banner) and `inst/stan/amm_eb_marginal_K.stan` and
  `inst/stan/amm_eb_conditional_K.stan` (Path B; bit-identical bodies to
  `amm_distrib_K.stan` modulo the `theta_ref_k_data` move + EB banner).
  The combined regime $K > 1 \wedge p > 1$ remains guarded with
  `gdpar_unsupported_feature_error` and is queued for Sub-phase 8.6.D.
* **(K, p)-aware dispatcher** in
  `.gdpar_eb_generate_stan_marginal()` /
  `.gdpar_eb_generate_stan_conditional()`: selects one of the four EB
  templates based on the resolved $(K, p)$ pair; the K = 1 + p = 1 leaf
  preserves the byte-identical 8.6.B template selection.
  `.gdpar_eb_check_stan_id_for_path(family, K, p)` enforces the
  per-path supported `stan_id` set: Path A supports
  $\{1, 2, 3, 4\}$ (inherits `amm_distrib_multi.stan`), Path B supports
  $\{1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13\}$ (inherits
  `amm_distrib_K.stan`'s full distributional-regression dispatch).
* **Three K-input patterns wired into `gdpar_eb()`**: the function now
  accepts the same three forms as `gdpar()` for declaring multi-slot
  models: (i) `formula` as `gdpar_formula_set` (e.g.,
  `gdpar_bf(y ~ a(x), phi ~ a(z))`); (ii) `amm` as a named list of
  `amm_spec`; (iii) `formula` with `a()`/`b()`/`W()` wrappers in its
  RHS. When the resolved K is > 1, the new internal orchestrator
  `.gdpar_eb_run_K()` runs Path B end-to-end.
* **Matrix-valued Proposition 7B\* correction** under $p > 1$
  (`.gdpar_eb_correction_matrix()`): returns
  $C^*_{g, \alpha} = \kappa(\alpha) \cdot \Sigma^{\text{marg}}_{\theta_{\text{ref}}}$
  as a $p \times p$ matrix; reduces algebraically to the scalar
  Proposition 7B of v07 Section 6 at $p = 1$. The `correction_applied`
  / `eb_correction_constant` slot of `gdpar_eb_fit` now carries either
  the scalar form (Sub-phase 8.6.B contract preserved bit-exactly) or
  the matrix form (Path A under D34).
* **`R/golden_helpers.R`**: `.gdpar_golden_roster_8_6_C()` canonizes
  17 configurations (10 Path A + 7 Path B; 13 fitable + 4 guard).
  `.gdpar_golden_eb_rds_name_8_6_C(family, path, regime)` produces the
  canonical filename
  `golden_eb_<family>_path<A|B>_<regime>.rds` so the 8.6.B and 8.6.C
  on-disk namespaces never collide.
* **`tests/testthat/test-eb_smoke_by_family.R`**: 17 new smokes appended
  (smokes 18-34). 4 Path A guard smokes run unconditionally (verify
  `gdpar_unsupported_feature_error` on Beta/Gamma/Student-t/Tweedie at
  $K = 1, p = 2$). 6 Path A fitable + 7 Path B fitable smokes are
  gated by `GDPAR_RUN_STAN_SMOKE_EB=1`.
* **`tests/testthat/test-eb_goldens_8_6_C.R`**: new structural test
  mirroring `test-eb_goldens_8_6_B.R`. Validates roster coverage,
  filename convention, Path A K=1+p>1 / Path B K>1+p=1 invariants,
  guard class canonicity, and manifest CSV cross-validation.
* **`data-raw/bootstrap_eb_goldens_C.R`**: new bootstrap script gated
  by `GDPAR_BOOTSTRAP_EB_GOLDENS=1` + `GDPAR_BOOTSTRAP_EB_GOLDENS_C=1`.
  Bootstraps the 13 fitable goldens with seed `20260525L` and 4
  metadata-only guard goldens; appends 17 rows to
  `tests/testthat/data/golden_manifest.csv` with `sub_phase = "8.6.C"`.
* **Decision D35 (emergent, no canonization reopened)**: three runtime
  corrections to `R/eb.R` necessary for the multivariate paths to run
  end-to-end. (i) Path-aware init dispatch in
  `.gdpar_eb_maximize_marginal()` (Path A / Path B delegate to
  cmdstanr's default unconstrained-space random init; K = 1 + p = 1
  preserves the explicit `.gdpar_eb_make_random_init()` bit-exactly).
  (ii) Path-aware extraction of `theta_ref` from Laplace draws:
  detects `theta_ref_k[j,k]` (Path B) and `theta_ref[j,k]` (Path A)
  variable name conventions. (iii) Reshape of `theta_ref_data` to a
  $J_{\text{groups}} \times p$ matrix in the K = 1 pipeline before
  the conditional HMC step (the Stan conditional template declares
  `array[J_groups] vector[p] theta_ref_data` and cmdstanr's
  auto-packing requires the R-side matrix shape under $p > 1$). All
  three are minimum-invasive patches at the R-Stan binding layer;
  none touches the canonized Stan templates or the Charter Section 2
  decisions.
* **Path B numerical sensitivity caveat**: families with strict
  $(0, 1)$ inverse-link (Beta) can exhibit conditional HMC instability
  under short warmup when the EB plug-in anchor is near the logit
  boundary. Workarounds: `iter_warmup >= 1000`, tighter `sigma_a_k`
  prior, or alternative anchor selection. Documented as deferred
  numerical caveat for Sub-phase 8.6.D Hessian-fragility expansion.
* Suite delta over Sub-phase 8.6.B baseline (2918 PASS): final
  **`[ FAIL 0 | WARN 0 | SKIP 82 | PASS 3021 ]`** (+103 PASS / +17
  SKIP), of which Sub-phase 8.6.C contributes +63 PASS / +17 SKIP on
  top of the +40 PASS canonized by 8.6.C pasos 1-4 (Session 11).

## Sub-phase 8.6.B — Empirical Bayes orchestrator base regime (2026-05-25)

* **`gdpar_eb(formula, family, amm, data, ..., eb_correction, laplace_control)`**:
  new exported function that fits an AMM canonical model via Empirical Bayes
  (Type II ML for `theta_ref` + conditional HMC for
  `xi = (a, b, W, sigma_*, phi)` given the plug-in estimate
  $\widehat\theta_{\mathrm{ref}}^{\,\mathrm{EB}}$). Sub-phase 8.6.B is the
  base regime $K = 1$, $p = 1$; multivariate $p > 1$ (8.6.C), multi-slot
  $K > 1$ (8.6.C), and mixture/hurdle families with `min_K > 1` (8.6.C) are
  rejected by the input guard with `gdpar_unsupported_feature_error`.
  Decisions 2.1, 2.2, 2.6, 2.7, 2.8 of `CHARTER_SUBFASE_8_6.md`.
* **Two dedicated Stan templates**: `inst/stan/amm_eb_marginal.stan` and
  `inst/stan/amm_eb_conditional.stan`. The marginal template is the
  cmdstanr::laplace() target of Step (i); the conditional template moves
  `theta_ref` from the parameters block to the data block (as
  `theta_ref_data`) and removes its prior so that Step (iii) HMC samples
  the conditional posterior of xi only. Both templates declare the
  dimension fields `K_slots` and `p_dim` so that Sub-phase 8.6.C can
  habilitate $K > 1$ and $p > 1$ without rewriting the data block.
  Decisions 2.3, 2.7 of `CHARTER_SUBFASE_8_6.md`.
* **Anti-fragility strategy of Charter §2.8 (four components)**:
  (i) preventive Hessian condition-number check
  (`laplace_control$kappa_threshold`, default `1e10`; aborts with
  `gdpar_eb_numerical_error` when violated);
  (ii) adaptive Levenberg-Marquardt ridge on the marginal covariance
  when its eigen-decomposition is singular or non-PD
  (`laplace_control$ridge_init`, default `1e-6`);
  (iii) multi-start optimize() with `laplace_control$multi_start_M`
  independent random inits (default `5`), retaining the init with the
  highest joint MAP log-density and reporting the dispersion across
  inits in `gdpar_eb_fit$diagnostics_numerical$multi_start_dispersion`;
  (iv) documented fallback that aborts with
  `gdpar_unsupported_feature_error` recommending FB when every init
  fails.
* **Proposition 7B scalar coverage correction** (default
  `eb_correction = TRUE`): the conditional credible intervals reported
  by `summary.gdpar_eb_fit` are inflated by a scalar factor derived
  from the marginal variance of `theta_ref`; setting
  `eb_correction = FALSE` issues a `gdpar_diagnostic_warning` advising
  of the expected $O(n^{-1})$ under-cover. Decision 2.6 of
  `CHARTER_SUBFASE_8_6.md`.
* **S3 methods**: `print.gdpar_eb_fit`, `summary.gdpar_eb_fit`,
  `print.summary.gdpar_eb_fit`, `coef.gdpar_eb_fit`,
  `predict.gdpar_eb_fit` (in-sample only in 8.6.B; `newdata` rejected
  with `gdpar_unsupported_feature_error` until 8.6.C).
* **Partial release of Sub-phase 8.6.B**: the smoke-by-family matrix
  of Charter §3.2 (17 cases), the 17 golden snapshots, the EB-vs-FB
  first-order coherence test, and the dedicated
  `test-eb_numerical_robustness.R` suite are listed as TODO in
  `HANDOFF_SUBFASE_8_6_B_PARCIAL.md` §5 and are deferred to the next
  clean session that closes Sub-phase 8.6.B.

## Sub-phase 8.5.B — External meta-learner comparator (2026-05-24)

* **`gdpar_compare_meta_learners(bridge, methods, newdata, data, seed)`**:
  new exported function that benchmarks a `gdpar_causal_bridge` against a
  user-supplied set of external meta-learner adapters on a common
  evaluation grid. Reports per-method `cate_mean` (and native CIs when the
  adapter exposes one) plus three concordance matrices (RMSE, Pearson,
  mean absolute discrepancy) over every ordered pair of methods.
  Descriptive by construction: no tests of hypothesis, no claims of
  algorithmic equivalence across methods of different inferential origin.
  Decisions 2.1, 2.3, 2.5, 2.6 of `CHARTER_SUBFASE_8_5_B.md`.
* **`gdpar_meta_learner_adapter(name, fit_predict_fun, predict_fun, ...)`**:
  constructor of the pluggable adapter contract. The mandatory
  `fit_predict_fun` carries the full fit-and-predict cycle on
  `(X, Y, T, X_newdata, level, seed_run)`; the optional `predict_fun`
  re-evaluates the meta-learner on a fresh grid by reusing the fitted
  `state` returned by `fit_predict_fun`, which lets
  `predict.gdpar_meta_learner_comparison` avoid a refit when supported.
  `is_gdpar_meta_learner_adapter()` and `print.gdpar_meta_learner_adapter()`
  round out the interface.
* **Two reference adapters distributed in the package**:
  `gdpar_adapter_grf()` wraps `grf::causal_forest` R-side with native CIs
  from `predict(..., estimate.variance = TRUE)` and the normal
  approximation; `gdpar_adapter_econml()` wraps
  `econml.dml.CausalForestDML` Python-side via `reticulate`, with native
  CIs from `effect_interval()`. Both expose `fit_predict_fun` and
  `predict_fun`; the Python reference inside the EconML adapter's `state`
  becomes invalid after `saveRDS` round-trip and the `predict_fun` aborts
  cleanly with `gdpar_unsupported_feature_error` in that case. Decision
  2.1 of `CHARTER_SUBFASE_8_5_B.md`.
* **S3 methods**: `print.gdpar_meta_learner_comparison` for a concise
  summary including the three concordance matrices;
  `summary.gdpar_meta_learner_comparison` for a structured object with
  long-format pairwise metrics, per-method ATE (with CI bounds derived
  from per-observation native CIs when available), and timing;
  `predict.gdpar_meta_learner_comparison` to re-evaluate the CATE on a
  new grid for every method, reusing cached state via the adapter's
  `predict_fun` and emitting a `gdpar_diagnostic_warning` when any
  adapter would require a full refit.
* **Scalar-outcome restriction**: bridges built from fits with `K > 1`
  or `p > 1` are rejected with `gdpar_unsupported_feature_error`;
  multi-output external adapters are queued for Block 9. Decision 3.8
  of `CHARTER_SUBFASE_8_5_B.md`.
* **New `Suggests`**: `grf` (R-side reference adapter) and `reticulate`
  (Python bridge for EconML). `econml` itself is a Python module the
  user installs in their active Python environment (e.g.
  `reticulate::py_install("econml")`); the package does not install
  Python dependencies automatically. Decision 4.6 of
  `CHARTER_SUBFASE_8_5_B.md`.
* **Two new vignettes**:
  `vignette("v08c_meta_learner_comparison")` canonizes the comparator
  contract, the concordance criterion, and the limits of cross-method
  comparison; `vignette("vop06_meta_learner_comparison")` is the
  step-by-step operational recipe (grf CRAN-valid; EconML with chunks
  gated by Python availability). Decision 2.2 of
  `CHARTER_SUBFASE_8_5_B.md`.

## Sub-phase 8.5.A — T-learner causal bridge, AMM-side (2026-05-24)

* **`gdpar_causal_bridge(fit_treat, fit_ctrl, newdata, type, level)`**: new
  exported function that estimates the conditional average treatment effect
  (CATE) from a pair of independent `gdpar_fit` objects fitted to the two
  arms of a treatment / control design. Implements the T-learner
  meta-learner of Kuenzel et al. (2019) on the AMM-side: each arm is fitted
  independently and the per-observation CATE is the difference of the two
  predictive distributions evaluated on a common evaluation set. Decision
  2.1 of `CHARTER_SUBFASE_8_5_A.md`.
* **Strict structural compatibility checks** on construction: the two fits
  must share the family identifier and link (or per-slot family identifiers
  when `K > 1`), `K`, `p`, AMM level, modulating basis type, anchor value,
  and the covariate column structure of every AMM component. Mismatches
  abort with `gdpar_unsupported_feature_error`. Hierarchical fits
  (`stan_data$use_groups == 1L`) are out of scope and likewise abort; the
  abort uses the same condition class as `gdpar_bvm_check` so user code
  that handles unsupported-feature errors covers both helpers.
* **S3 methods**: `print.gdpar_causal_bridge` for a concise object summary;
  `summary.gdpar_causal_bridge` for a per-observation CATE table plus the
  marginal ATE with credible bounds; `predict.gdpar_causal_bridge` for
  re-evaluating the CATE on a fresh `newdata` without re-running the
  compatibility checks.
* **New canonical vignette** `vignette("v08b_cate_ite_bridge_implementation")`
  canonizing the T-learner AMM-side: definition, identification assumptions
  inherited from `v02` plus residual no-confounding, estimator,
  per-observation credible bounds by posterior quantiles, identifiability
  per arm, a minimal reproducible example, limitations of the T-learner,
  and open questions (O*-CATE) for Block 9. Continues the CATE/ITE
  positioning of `vignette("v08_cate_ite_positioning")`.
* **No new heavy `Suggests`**: the comparator against external meta-learners
  (`grf`, `causalForest`, EconML via `reticulate`) is queued for Sub-phase
  8.5.B and adds no dependencies here.

## Sub-phase 8.3.10 — release-gate of Session 8.3 (2026-05-22)

* **`coef.gdpar_fit` for K-individual fits** (`K > 1`, distributional
  regression on `amm_distrib_K.stan`). Decision E4.A of sub-phase
  8.3.10: returns a named list of length `K` whose entries are
  `gdpar_coef` objects (each with `p = 1L`), one per slot. The modulating
  block `W_raw` is globally shared across slots in the K-individual
  template (canonical decision "Scope of W: global", handoff 28); the
  resulting `W` component is replicated identically across every slot's
  `gdpar_coef` when any slot declared a non-NULL `W`. K-individual
  fits with grouping (`J_groups > 1`) raise
  `gdpar_unsupported_feature_error` and remain queued for Session 8.4.
* **`predict.gdpar_fit` with `newdata` for K-individual fits**. Mirrors
  `predict_from_newdata_multi` per slot via the new internal helper
  `predict_from_newdata_K`. Each slot's `eta_k` is rebuilt from
  `theta_ref_k[k]`, `a_coef_k[k]`, `c_b_k[k]`, and the globally-shared
  `W_raw` + `sigma_W` (with the basis-difference vector evaluated at
  the slot-specific `theta_ref_k[k]` and `theta_anchor_K[k]`). The
  polynomial `W` branch is supported; B-spline `W` on new data remains
  queued for Session 8.4 (the in-sample path through
  `apply_W_basis_diff()` in Stan continues to support both). Grouping
  (`J_groups > 1`) on new data is also queued for Session 8.4.
* **Two new vignettes** added for pedagogical coverage of the
  intermediate AMM regime (decision E5.C of sub-phase 8.3.10):
    - `vignette("vop04_amm_intermediate")` — B-spline `W` bases and
      heterogeneous families per slot.
    - `vignette("vop05_distributional_K_dharma")` — distributional
      regression `K > 1` API, residual diagnostics G1 / G2 / G3, and
      optional DHARMa integration via `gdpar_dharma_object()`.
* **Bit-exact compare-path tests for the 14 8.3.9 goldens**, tiered by
  computational cost (decision E6.C of sub-phase 8.3.10):
    - Tier 1 (`GDPAR_GOLDEN_8_3_9_COMPARE_TIER1=1`, ~9 min): 12 configs
      across sub-phases 8.3.4 / 8.3.5a / 8.3.6 / 8.3.7 / 8.3.8 / 8.3.9.
    - Tier 2 (`GDPAR_GOLDEN_8_3_9_COMPARE_TIER2=1`, ~40 min): the
      Tweedie K=3 config (~38 min by itself) and `bspline_p2`.
    - The bootstrap script `inst/extdata/scripts/bootstrap_8_3_9_goldens.R`
      is now `source()`-able from the test file; its top-level driver
      runs only when invoked as a top-level Rscript with
      `GDPAR_BOOTSTRAP_8_3_9=1`.
* **Test suite at the close of 8.3.10**: 2630 PASS / 0 FAIL / 0 WARN /
  40 SKIP (delta vs 8.3.9 baseline: +55 = 32 coef K>1 + 23 predict
  newdata K>1; SKIP unchanged).
* **Documentation regenerated** (`devtools::document()`): `coef.gdpar_fit.Rd`
  updated to describe the named-list return for K > 1; no NAMESPACE
  changes.

## Sub-phase 8.3.9 — residual diagnostics G1 / G2 / G3 and golden bootstrap (2026-05-22)

* **Three residual layers** (decision D4 of sub-phase 8.3.9):
    - **G1** — deviance and Pearson residuals per family (closed-form
      where available; Pearson-like fallback for mixtures and Hurdle).
    - **G2** — Bayesian randomized quantile residuals (Dunn-Smyth 1996)
      via the empirical CDF of posterior-predictive draws, jittered
      for discrete responses and mapped through `qnorm()`.
    - **G3** — posterior-predictive checks (PPC) via the new S3 method
      `pp_check.gdpar_fit()` off the `bayesplot::pp_check` generic
      (five PPC types: `dens_overlay`, `hist`, `ecdf_overlay`, `stat`,
      `intervals`).
* **`residuals.gdpar_fit()`** exported with `type ∈ {"quantile",
  "response", "pearson", "deviance"}`; `gdpar_posterior_predict()`
  exported for raw posterior-predictive draws; `gdpar_dharma_object()`
  exported for optional DHARMa integration (DHARMa is a `Suggests`
  dependency, decision E1.A; the package's minimalist `Imports` is
  preserved).
* **14 golden RDS snapshots** bootstrapped under cmdstan 2.38.0
  (decision D1=(1D)): lognormal `K=2`, Student-t `K=3`, Tweedie `K=3`,
  ZIP / ZINB / Hurdle-Poisson / Hurdle-NB, heterogeneous Gauss+Beta /
  Gauss+Gamma / NB+Beta (`K=2`), Gaussian `K=1` bspline `p=1` / bspline
  `p=2` / polynomial `p=2`, and a heterogeneous Gauss+Beta `K=2` bspline
  config. Manifest CSV with 16 columns of reproducible metadata
  (SHA256 `fit_code_hash`, seed, sub-phase, cmdstan / R / DHARMa
  versions, etc.) at `tests/testthat/data/golden_manifest.csv`
  (decision E3.A).

## Sub-phase 8.3.8 — B-spline `W` bases (2026-05-22)

* **`W_basis(type = "bspline", knots, degree, boundary_knots)`** added.
  Stan-side Cox-de Boor recursion via `bspline_basis_eval()` and the
  helper `apply_W_basis_diff()`, used uniformly across the three Stan
  templates (`amm_main.stan`, `amm_distrib_multi.stan`,
  `amm_distrib_K.stan`).
* **Backward-compat preserved**: polynomial `W` paths are
  bit-exact-preserved by constant propagation; the data field
  `W_is_polynomial` is replaced by `W_type_id` + `W_n_knots_full` +
  `W_knots_full` + `W_degree`.

## Sub-phase 8.3.7 — heterogeneous families per slot (2026-05-21)

* **Named-list `family` argument for `gdpar()` with `K = 2`**
  (decision D5). Slot 1 is the canonical primary location family
  (`gaussian`, `poisson`, `neg_binomial_2`, `beta`, `gamma`,
  `lognormal_loc_scale`); slot 2 may declare any compatible auxiliary
  family with coherent `support`.
* **`apply_inv_link_by_id()` Stan helper** added; refactors the K=2
  branches (stan_id 1 / 3 / 5 / 6 / 7) to dispatch the per-slot
  inverse link uniformly. Goldens K=2 Beta + Gamma were re-bootstrapped
  to absorb the refactor (decision D7=(a)).

## Sub-phase 8.3.6 — mixture likelihoods (2026-05-21)

* **ZIP `K=2`** (stan_id 10, `mu log + pi logit`),
  **ZINB `K=3`** (stan_id 11, `mu log + phi log + pi logit`),
  **Hurdle-Poisson `K=2`** (stan_id 12), and
  **Hurdle-NB `K=3`** (stan_id 13). Hurdle truncation enforced via
  `log1m_exp(neg_binomial_2_log_lpmf(0 | mu, phi))` for numerical
  stability.

## Sub-phase 8.3.5 — `K = 3` likelihoods (2026-05-21)

* **Student-t `K=3`** (stan_id 8, `mu identity + sigma log + nu log`).
* **Tweedie `K=3`** (stan_id 9, `mu log + phi log + p identity bounded
  (1.01, 1.99)`). Custom `tweedie_lpdf` in the `functions{}` block via
  a hybrid Dunn-Smyth series + Wood (2017) saddlepoint switchover at
  `tau = 0.4`; `tweedie_rng` via compound Poisson-gamma.

## Sub-phase 8.3.4 — bi-parametric `K = 2` likelihoods (2026-05-21)

* **Beta `K=2`** (stan_id 5, `mu logit + phi log`),
  **Gamma `K=2`** (stan_id 6, `mu log + shape log`),
  **log-normal `K=2`** (stan_id 7, custom-family wired via descriptor).

## Sub-phase 8.3.3 — brms-style formula API for `K > 1` (2026-05-20)

* **`gdpar_bf(y ~ a(x1), sigma ~ a(x2))`** constructor + **`gdpar_formula_set`** S3 class. Three equivalent input forms canonicalise to the same internal representation.
* **`individual_params` argument removed** from `gdpar_family()` (decision 4C superseded by canonisation D). `gdpar_family()` now emits *all* eligible `param_specs`; the input to `gdpar()` is the single source of truth for `K`.

## Breaking change in the AMM canonical implementation (2026-05-18)

* **Centering conditions (C2) and (C3) of Block 1 are now enforced
  empirically by column-wise centering of the design matrices Z_a and
  Z_b in the R-side AMM design constructor (`amm_spec_designs()`,
  `amm_spec_multi_designs()`), and no longer by a sum-to-zero
  reparametrization on the basis coefficients inside Stan.**
* Rationale: with `colMeans(Z_a) = 0` enforced empirically, the
  expectation `E_mu[a(X)] = colMeans(Z_a) * a_coef` vanishes for any
  `a_coef`, so (C2) is satisfied for the full `J_a`-dimensional space
  of basis coefficients. The previous implementation enforced the
  same centering twice (in R and again in Stan), and the second
  enforcement, sum-to-zero on `a_coef`, was a strictly stronger
  restriction not derived from (C2): it confined `a_coef` to the
  `J_a - 1`-dimensional subspace `{a : sum(a) = 0}` for every fit.
* Empirical consequence: data-generating truths with effects that do
  not have zero coefficient sum (e.g. `a_coef = (1.0, -0.8)` on two
  continuous covariates) are now recovered correctly. Previously they
  were projected onto the sum-zero subspace at fit time, inflating
  residual variance and degrading ELPD. The S7 (`heterog_p3`) scenario
  of the Block 6 adversarial bench, where gdpar previously drifted
  88 ELPD points outside `2*se` versus the other four competitors, now
  cleanly enters the indistinguishable cluster (gdpar = -1378.45 vs
  best mgcv = -1375.50, delta 2.95 within `2*se = 70`).
* Affected files: `inst/stan/amm_main.stan`,
  `inst/stan/amm_main_multi.stan`, `R/stan_codegen.R`
  (`generate_a_blocks_multi` NCP/CP/mixed variants), `R/preflight.R`,
  `R/preflight_multi.R`. Documentation revised in
  `vignette("v01_amm_identifiability")` Section 10.1 and
  `vignette("v03_special_cases")` Sections 6 and 8. Roxygen for
  `amm_spec()` and `gdpar_prior()` reworded. The Block 6 multi golden
  fits (`smoke_p2_auto.rds`, `smoke_p2_cp.rds`, `smoke_p2_ncp.rds`)
  were regenerated against the new codegen.
* Test suite: 1099 PASS / 0 FAIL / 0 WARN / 12 SKIP without
  regressions (one additional `expect_false` for the absent
  sum-to-zero emission on the scalar template).

## API change (2026-05-19)

* `gdpar()` gains argument `id_check_rigor = c("full", "fast")` with
  default `"full"` (backward-compatible). The flag is propagated to
  `.gdpar_multi()` and to `gdpar_check_identifiability(..., rigor =
  id_check_rigor)`. The `"fast"` value downgrades the C4-bis
  cross-coordinate check to a single consolidated warning at the end
  of the per-coordinate loop instead of aborting, enabling fits with
  intentional overlap between the additive and modulating channels
  combined with informative priors. Documented in `?gdpar` and in
  `vignette("v01_amm_identifiability", package = "gdpar")` Section
  6.6.1.7.

## Block 6 — synthetic adversarial validation TOP4 with coord-wise factorization (in progress, 2026-05-18 / 2026-05-19)

### New features

* `gdpar_loo()` exported (flagged `@keywords experimental`): PSIS-LOO
  approximate cross-validation via `loo::loo()` on the per-observation
  log-likelihood emitted by the Stan template. Two aggregations:
  `"subject"` (default; sums the coord-wise log-likelihoods to
  per-subject scale, matching `brms::brm()` multivariate fits with
  `set_rescor(FALSE)`) and `"cell"` (treats each `(i, k)` pair as an
  independent observation, useful for per-coordinate Pareto-k
  diagnostics).
* Multi adversarial bench harness under `inst/benchmarks/`:
  * `inst/benchmarks/synthetic_hard_multi.R`: 9 scenario generators
    (S0..S8) covering sanity, adversarial, and stress tiers with
    `n_train = 800`, `n_test = 200` per scenario. Includes a
    Bernoulli `p = 2` adversarial scenario (S5) and a Gaussian
    `p = 2` residual-correlation scenario (S8) where `brms` is run a
    second time with `set_rescor(TRUE)` under the label `brms_rescor`
    to verify recovery of the true residual correlation.
  * `inst/benchmarks/scripts/bench_competitor_<method>.R` wrappers
    for gdpar, brms, mgcv, INLA, and rstanarm. PSIS-LOO is computed
    on principled posterior draws for all five methods: gdpar /
    brms / rstanarm use their native Stan posteriors; mgcv uses the
    Bayesian interpretation of REML (Wood 2017, sec 6.10) via
    `mvnfast::rmvn(mu = coef(gam_obj), sigma = vcov(gam_obj))` to
    sample from the gaussian Laplace approximation; INLA uses
    `inla.posterior.sample()` with a CPO-sum cross-check that flags
    `inla_cpo_loo_divergent` when the two values disagree by more
    than `2 * se_elpd`.
  * `inst/benchmarks/run_synthetic_hard_multi.R` runner, env-gated
    by `GDPAR_BENCH_MULTI=1`, writes the long-tidy CSV
    `synthetic_hard_multi_results.csv` with columns including
    `max_pareto_k`, `n_pareto_k_above_07`, `convergence_flag`, and
    `recovered_rescor` (NA except for `brms_rescor` rows in S8).
  * `inst/benchmarks/scripts/report_synthetic_hard_multi.R`
    visualizer producing `synthetic_hard_multi_plot.png`,
    `synthetic_hard_multi_wall.png`,
    `synthetic_hard_multi_pareto.png`, and a markdown verdict table
    `synthetic_hard_multi_summary.md`.
* `vignettes/vop02_arbitrary_p.Rmd` extended with Section 10
  "PSIS-LOO via `gdpar_loo()`" documenting both aggregations,
  Pareto-k caveats, and the experimental-status flag.

### DESCRIPTION

* `Suggests:` adds `brms`, `INLA`, `rstanarm`, `scoringRules`,
  `gratia`, `mvnfast` for the Block 6 bench harness. All entries are
  optional and gated by `requireNamespace()` at call sites.
* `Additional_repositories:` adds
  `https://inla.r-inla-download.org/R/stable` so that `INLA` can be
  installed from CRAN-compatible tooling.

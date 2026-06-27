# Block 9.2.S — Synthetic recovery benchmark (known ground truth)

**Status: COMPLETE — coverage gate satisfied (3200/3200 cells by status).**
Generated B9.18 (2026-06-01). Roster: gdpar (subject) + mgcv + brms + INLA +
rstanarm. Driver `run_recovery_9_2_S.R`; scenarios `scenarios_recovery_factorial.R`;
adapters `scripts/recovery_adapters.R`; metrics `recovery_metrics.R`. Per-scenario
matrix in `matrix_64x5_rmse_mean.csv`; per-cell manifest `manifest_9_2_S.csv`
(3200 rows); cell objects `cells_9_2_S/*.rds`.

No result here is publishable from partial coverage: it is reported only now
that every one of the 3200 cells carries a recorded status.

---

## 1. Design

A `2 x 2^5 = 64`-scenario factorial with **a-priori known structural truth**, so
each method's estimated conditional moments are scored against the true mean,
scale and 0.95 tail quantile on a fixed evaluation grid (n_grid = 200). Two
response regimes, each a 2^5 factorial over five on/off pathologies:

| bit | regime C (continuous) | regime K (count) |
|-----|-----------------------|------------------|
| b1  | multicollinearity     | multicollinearity |
| b2  | nonlinearity          | nonlinearity      |
| b3  | autocorrelation       | autocorrelation   |
| b4  | heteroscedasticity    | zero-inflation    |
| b5  | heavy tails           | overdispersion    |

b1/b2/b3 are covariate-side and shared; b4/b5 are response-side and
regime-specific. The active response bits fix the gdpar family/K:
gaussian K1 / gaussian K2 (heterosc.) / student_t K3 (b5) for regime C; poisson
K1 / neg_binomial_2 K2 (b5) / zip K2 (b4) / zinb K3 (b4+b5) for regime K.

- **R = 10** repetitions per scenario; **n = 2000** (1600 train / 400 test);
  truth seeds 91005..91068 (deterministic, derived per-cell, not frozen-golden).
- **5 methods x 64 scenarios x 10 reps = 3200 cells.** HMC budget for the
  Bayesian methods: 2 chains, 500 warmup, 1000 sampling.
- Metrics (criterion 2.S, DESIGN 2.5): `rmse_mean` and `coverage90` are global;
  `rmse_sd` grades scale recovery (heteroscedasticity / overdispersion);
  `tail_metric` grades the 0.95 quantile (heavy tails); on **autocorrelation the
  only requirement is finiteness** (robustness, not leadership) because the true
  mean is AR-free and gdpar does not model AR. `elpd_loo` is a secondary
  cross-check only.

### 1.1 Nonlinearity: real comparison (Philosophy F, D82 + D83)

Each method recovers `sin(x1)` with **its own native smoother**, on equal footing
— so b2 is a *real comparison*, not a misspecification stratum:

- **mgcv / brms** — their native penalized `s(x1)`.
- **gdpar** — `a(bs(x1))` in the mean (its native `model.matrix` path).
- **INLA** — its native second-order random walk `f(x1, model="rw2")` (**D83**;
  see below).
- **rstanarm** — the only method with no native penalized smoother — consumes a
  **frozen, shared B-spline basis** of x1 (df = 9, degree 3, fit on train and
  applied identically to the grid) as plain fixed effects.

**D83 (B9.18).** The original D82 gave INLA the same frozen fixed-effect basis as
rstanarm. That broke INLA's Newton-Raphson optimizer for the Student-t family
(heavy tails + collinear flexible mean -> "lambda < 1/lambda_lim"; regularizing
the fixed-effect prior did not help). INLA *does* have a native penalized
smoother (rw2); bucketing it with rstanarm was a category error, the same kind
corrected in B9.17 for gdpar's W(theta_ref). Switching INLA to rw2 converges
across all families. A robustness fallback (Gaussian latent + empirical-Bayes
hyperparameters) fires only when the default INLA strategy diverges (1 cell in
3200: the hardest count realization, zip + nonlinearity + AR + multicollinearity).

### 1.2 The gdpar-DGP-aware asymmetry (honest)

gdpar's K>=2 slot specifications are **DGP-aware**: each slot models the
dependence the truth actually has (sigma~a(z) under heteroscedasticity; pi~a(z)
under zero-inflation; nu, phi intercept-only since the truth never modulates
them with covariates). Competitors run at their **defaults**. Where a competitor
can express the same structure (mgcv `gaulss`/`scat`/`ziP`, brms distributional)
it does; where it cannot (INLA/rstanarm have no distributional scale; rstanarm
no Student-t / zero-inflation) it is recorded `unsupported`. This asymmetry is
the whole point of the response-side strata: it measures what explicit
distributional modelling buys.

---

## 2. Coverage gate — 100% by status

| method   | ok  | unsupported | genuine capability gap |
|----------|-----|-------------|------------------------|
| gdpar    | 640 | 0           | — (fits the whole grid) |
| brms     | 640 | 0           | — |
| mgcv     | 560 | 80          | zinb (no zero-inflated NB) |
| inla     | 560 | 80          | heteroscedastic gaussian (no distributional scale) |
| rstanarm | 240 | 400         | student_t, zip, zinb, heteroscedastic gaussian |
| **total**| **2640** | **560** | **3200 / 3200 = 100%** |

All 560 `unsupported` are genuine package gaps (robustness passes, DESIGN 2.5.4);
**zero `fit_failed`, zero `metrics_failed`.**

---

## 3. Robust reading: median, plus two catastrophic failures reported separately

The cross-cell distribution of `rmse_mean` is heavy-tailed: **exactly two cells**
(of 3200) blow up, both fully-Bayesian GLM-style competitors on a single count
realization where the posterior predictive draws explode:

| method   | cell                    | family         | rmse_mean |
|----------|-------------------------|----------------|-----------|
| brms     | K_10010 rep05           | neg_binomial_2 | 2.4e10    |
| rstanarm | K_11100 rep10           | poisson        | 1.4e4     |

**gdpar, mgcv, INLA: zero explosions** (gdpar max rmse_mean = 0.41; mgcv 0.48).
The mean is therefore not a usable central tendency; all leaderboards below use
the **median** and the two explosions are reported as a robustness finding, not
folded into the score (three-layer rigor: robust central tendency separate from
the tail-failure count).

---

## 4. Global leaderboard (median over each method's ok cells)

> Caveat: each method's median is over a **different cell subset** (gdpar/brms 64
> scenarios; mgcv/INLA 56; rstanarm 24). Cross-method global medians are
> indicative only. The rigorous comparisons are the per-scenario winner tally
> (Sec. 5) and the same-cell strata (Sec. 6).

| method   | rmse_mean | rmse_sd | tail_metric | coverage90 | wall (s, median) |
|----------|-----------|---------|-------------|------------|-------------------|
| gdpar    | 0.0599    | 0.1043  | 0.3464      | 0.850      | 171.4 |
| mgcv     | 0.0597    | 0.1320  | 0.3937      | 0.878      | 0.1 |
| brms     | 0.0720    | 0.1327  | 0.4637      | 0.820      | 96.0 |
| inla     | 0.1068    | 0.4688  | 0.9579      | 0.685      | 0.5 |
| rstanarm | 0.0906    | 0.0783* | 0.3953      | 0.885      | 6.4 |

\* rstanarm's low rmse_sd is over its 24 simple ok scenarios only (no
heteroscedastic / heavy-tailed cells), not comparable to the full-grid methods.

---

## 5. Per-scenario winner tally (lowest median rmse_mean)

Counting, for each of the 64 scenarios, which method recovers the mean best
**among the methods that fit that scenario**:

- **Full roster {gdpar, mgcv, brms}** (the three that fit all 64): **gdpar 36 /
  brms 16 / mgcv 12.** By regime — continuous: gdpar 15, brms 11, mgcv 6;
  **count: gdpar 21, brms 5, mgcv 6** (gdpar wins 21/32 of the count grid).
- **All five methods** (best among whoever fits): gdpar 33, brms 12, mgcv 9,
  rstanarm 6, inla 4.

gdpar wins the most scenarios overall and dominates the count regime.

---

## 6. Strata (criterion 2.S) — same-cell comparisons

### 6.1 gdpar's declared turf (response-side distributional structure)

| stratum | metric | gdpar | mgcv | brms | inla | rstanarm |
|---------|--------|-------|------|------|------|----------|
| heavy tails (student_t, C b5=1) | tail_metric | **0.228** | 0.621 | 0.620 | 1.153 | n/s |
| zero-inflation (count b4=1) | rmse_mean | **0.118** | 0.262 | 0.256 | 0.855 | n/s |
| heteroscedasticity (gaussian K2, C b4=1 b5=0) | rmse_sd | **0.064** | 0.067 | 0.064 | n/s | n/s |
| overdispersion (count b5=1) | rmse_sd | 0.198 | **0.186** | 0.224 | 0.575 | 0.214 |

(n/s = unsupported.) gdpar **dominates tail-quantile recovery** (2.7x better than
the next method) and **zero-inflated count means** (2.2x better) — exactly the
structure its explicit distributional slots model. It is co-best on
heteroscedastic scale (with brms; INLA/rstanarm cannot fit it at all) and a close
second to mgcv's `nb` on overdispersion SD.

### 6.2 Robustness strata (NOT gdpar's claimed turf)

| stratum | metric | gdpar | mgcv | brms | inla | rstanarm |
|---------|--------|-------|------|------|------|----------|
| autocorrelation (b3=1) | rmse_mean | 0.084 | 0.083 | 0.113 | 0.147 | 0.144 |
| nonlinearity (b2=1) | rmse_mean | 0.088 | **0.077** | 0.098 | 0.137 | 0.129 |

- **Autocorrelation** — gdpar does not model AR; it stays finite, well-controlled
  and is top-2 (essentially tied with mgcv). Requirement (robustness, not
  leadership) is met. No method recovers the AR-inflated uncertainty well
  (coverage90: gdpar 0.75, mgcv 0.71, others 0.53-0.65 — gdpar least bad).
- **Nonlinearity** — under the real comparison (D82/D83), gdpar **loses narrowly
  to mgcv's penalized spline** (0.088 vs 0.077) and beats the rest. This is the
  honest, anticipated outcome: a flexible mean in x1 is not gdpar's claimed
  advantage, and a dedicated penalized-spline GAM edges it. gdpar can lose b2.

### 6.3 Calibration (median coverage90 by stratum)

gdpar's 90% mean-interval coverage is stable at 0.84-0.85 across strata (slightly
conservative-to-nominal); mgcv similar; INLA badly undercovers the count
zero-inflation cells (0.05). gdpar is the most consistently calibrated of the
fully-Bayesian methods.

---

## 7. Compute cost (honest trade-off)

Median wall time per fit: **gdpar 171 s**, brms 96 s, rstanarm 6.4 s, INLA 0.5 s,
**mgcv 0.1 s**. gdpar is by far the most expensive (consistent with Block 7:
gdpar buys its distributional accuracy at ~1000x mgcv's cost). For a point
mean-recovery on a simple cell, mgcv is overwhelmingly cheaper for an
indistinguishable result; gdpar's cost is justified only where its
distributional modelling wins (Sec. 6.1).

---

## 8. Honest verdict (criterion 2.S, no overreach)

- gdpar **leads where it is designed to**: heavy-tail quantiles and zero-inflated
  count means (decisive), heteroscedastic scale (co-best), overdispersion (close
  second). It wins the most scenarios overall (33/64) and dominates the count
  regime (21/32 full-roster).
- gdpar is **robust where it makes no claim**: finite and top-2 under
  autocorrelation it does not model; zero catastrophic failures across 640 fits
  (brms and rstanarm each explode once).
- gdpar **can lose, honestly**: it is a narrow second to mgcv on pure
  nonlinearity (not its turf), and it is the slowest method by a wide margin.
- The win on the response-side strata reflects the **DGP-aware asymmetry** of
  Sec. 1.2: gdpar is given slots matched to the truth, the competitors that can
  express the same structure are too, and the ones that cannot are recorded
  `unsupported`. The benchmark measures what explicit distributional modelling
  buys — not a claim that gdpar is uniquely capable.

gdpar is competitive-to-dominant on distributional recovery, robust elsewhere,
and expensive. That is the recovery story 9.2.S supports — no more, no less.

---

## 9. Knobs / caveats (revisable)

- `.rec_smooth_df = 9` (frozen B-spline for gdpar/rstanarm) ~ mgcv k=10.
- `.rec_inla_rw2_groups = 50` (INLA rw2 discretization; effective df chosen by
  INLA's hyperprior, D83).
- `elpd_loo` recorded for gdpar (480/640; `gdpar_loo` returns NULL on some K>=2
  count cells), brms (640), rstanarm (240); skipped by design for mgcv/INLA
  (DESIGN 2.4). Secondary cross-check only; not used in the leaderboards.
- Global per-method medians (Sec. 4) are over different cell subsets; use the
  per-scenario tally and same-cell strata for method comparison.

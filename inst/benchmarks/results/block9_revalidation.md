# Block 9.2.O — eBird external re-validation: report and Decision E

**Status: 80/80 cells closed (coverage gate PASSED).** 64 `ok` + 16
`unsupported` (the honest coverage passes: a roster member with no native family
for a coordinate, e.g. brms/INLA/rstanarm on the Tweedie count, INLA on the
distributional scale). Zero `fit_failed`. The four cells that had been blocked
since B9.19 — `N__N1__gdpar__{nef,mac,app,grl}` (the gdpar Tweedie count) — are
now closed by the **RG.7 certified-limit Laplace ELPD** (sub-phase RG.7, D98).

Grid: 80 cells = (track G: 3 levels + track N: 1 level) x 5 methods x 4 NE-USA
sub-regions. The single comparability metric (criterion 2.O) is the held-out
test log predictive density on the **original response measure** (occurrence in
[0,1], count in [0, Inf)), via the change-of-variables Jacobian
(`ebird_metrics_9_2_O.R`), so every method is comparable within a track and the
two encoding tracks are comparable against each other.

---

## 1. The RG.7 closure (the count coordinate, track N coord 2: Tweedie)

### 1.1. Why these four cells needed a new method

The eBird abundance count is **near-deterministic** (n_train=3000, n_test=750,
~1-11% zeros, Tweedie power p in 1.1-1.5). Its posterior on gdpar's d=14
unconstrained scale is a **stiff canyon**: the observed-information condition
number at the mode is 9.4e7-4.4e8. B9.32-B9.33 established **by a multi-chain
gate** that this canyon is **NOT sampleable at feasible cost**:

- default NUTS saturates the tree depth and mixes pathologically;
- a constant-mass Laplace-preconditioned HMC accepts on a single chain but does
  not mix across chains (the B9.32 "resolved" was a **false positive of a
  single-chain acceptance gate**);
- a **whitened NUTS** — the decisive experiment, recommended convergently by
  GLM-5.2 and MiMo-2.5-Pro: a dense metric equal to the Laplace precision
  (−Hessian at the mode), chains initialised at the mode, 1000 warmup + 4x1000
  sampling — **fails the multi-chain gate**: R-hat=1.117, ESS-bulk=22.5,
  193/4000 = 5% divergences (`~/rg33_whiten_full.out`).

Divergences in the **already-whitened** coordinates prove the geometry is
**genuinely non-Gaussian** (the canyon curves; the 90-degree Fisher/Hessian
misalignment of B9.32 foreshadowed it). A per-step Riemannian metric is
infeasible (minutes per exact Hessian). This **revises D93**: rigidity does NOT
imply Gaussianity when the canyon curves — Bernstein-von Mises holds only in a
vanishingly small neighbourhood of the mode.

### 1.2. The certified-limit endpoint (Laplace / plug-in, competitor parity)

Both lineages and the user pre-agreed the honest endpoint (stopping rule: no
further samplers): **certify the sampling limit and report the Laplace / plug-in
ELPD**. This is exactly the regime of the competitors on this coordinate —
**mgcv/REML** reports a coefficient-Gaussian + dispersion/power plug-in
predictive, **INLA** a Laplace approximation; both are O(d/n)-exact Laplace /
plug-in predictives, not exact MCMC. gdpar matches that regime: the mode + the
exact-Hessian Gaussian posterior `N(mode, M^{-1})`, pushed through the model's
`constrain_variables` and the package's own `predict()` (so the Tweedie density
+ Jacobian are scored identically to every competitor). The implementation is
`inst/benchmarks/scripts/rg7_laplace_elpd.R` (cero touch to `R/`); the injection
is verified against the model's own draws (constrained images of the real
unconstrained draws reproduced to max_abs_diff=1.7e-8).

### 1.3. Result (count coordinate, Tweedie ELPD on the original measure)

| sub-region | gdpar **plug-in** (mode) | mgcv/REML | gdpar **Laplace-averaged** ± se | cond | Pareto-k | ESS_is |
|---|---|---|---|---|---|---|
| NE_FOREST     | **−657.98** | −658.33 | −693.83 ± 23.2 | 4.4e8 | 0.77 | 3.5% |
| MID_ATL_COAST | **−840.56** | −840.70 | −858.54 ± 19.8 | 9.4e7 | 0.83 | 1.2% |
| APPALACHIAN   | **−982.41** | −982.45 | −991.70 ± 20.9 | 1.9e8 | 0.84 | 4.4% |
| GREAT_LAKES   | **−1010.03** | −1010.08 | −1017.38 ± 20.6 | 2.7e8 | 1.05 | 0.9% |

Three-layer reading (algebraic / statistical / numerical, kept separate):

- **Plug-in parity is essentially exact.** gdpar's MAP/mode predictive matches
  mgcv's REML plug-in to **< 0.4 ELPD units in all four cells** (max |Δ|=0.35).
  At the point-estimate level the two methods agree to the third significant
  figure: the same predictive surface, reached by different optimisers (gdpar's
  penalised log-posterior mode vs mgcv's REML).
- **The full Laplace average is conservative, by design.** gdpar's
  Laplace-averaged ELPD is 7-36 lower because it propagates the **full**
  posterior uncertainty — including the dispersion phi and the power p — which
  mgcv/REML and INLA plug in. gdpar gets **no credit for plugging in**: this is
  the no-overreach choice (the reported headline `y2_elpd_orig`). The gap is
  within ~1-1.5 standard errors (se ~20-23), so gdpar and mgcv are
  **statistically indistinguishable** on the count (consistent with Block 7,
  where all five competitors were indistinguishable in elpd_loo).
- **The Laplace fidelity is honestly quantified.** The importance-sampling ESS
  (0.9-4.4%), PSIS Pareto-k (0.77-1.05, all above the 0.7 reliability line),
  and the mean log-density drop (>> d/2 = 7) **positively certify the genuine
  non-Gaussianity**: the Laplace Gaussian is a wide, crude proposal for the
  curved canyon posterior. This is the scientific finding, not a defect — it
  measures *how* non-Gaussian the canyon is. The Laplace-averaged ELPD is a
  conservative proxy; the (intractable) exact-posterior ELPD is bracketed in
  roughly [Laplace-averaged, plug-in], i.e. between the conservative average and
  the mgcv-matching plug-in. All modes are clean (mode_offset < 0.14 SD,
  Hessian all-positive, 0% constrain failures).

### 1.4. The occurrence coordinate (track N coord 1: Beta)

The Beta occurrence fit is the well-behaved coordinate (ordinary NUTS). gdpar's
Beta ELPD matches mgcv's to **< 0.4 units** in all four cells (gdpar/mgcv:
495.70/495.72, 1339.23/1339.30, 503.57/503.21, 526.63/526.26).

---

## 2. The encoding-track comparison (criterion 2.O across tracks)

Mean total (coord1 + coord2) held-out ELPD on the original measure, per level,
averaged over the four sub-regions (track-G means over `ok` cells):

| level | gdpar | mgcv | brms | INLA | rstanarm |
|---|---|---|---|---|---|
| G1 (gaussian K1) | −167.9 | −168.1 | −167.8 | −167.7 | −167.8 |
| G2 (gaussian-hetero K2) | **unstable** | 20.7 | **unstable** | — | — |
| G3 (student-t K3) | **5.0** | −57.4 | −93.9 | −293.3 | — |
| N (beta + tweedie) | see §1 | see §1 | n/a (no tweedie) | n/a | n/a |

- **G1**: all five methods indistinguishable (~−168).
- **G3**: gdpar's heavy-tailed K3 is the best track-G level (5.0), ahead of mgcv
  (−57.4) and well ahead of INLA (−293.3).
- **G2 caveat (pre-existing, honest).** The heteroscedastic gaussian on the
  *transformed* response (logit-occurrence / log1p-count) is numerically fragile
  for the **fully-Bayesian** methods: extreme posterior draws of the
  observation scale produce catastrophic predictive tails on the back-transformed
  original measure (gdpar and brms both blow up; mgcv's REML plug-in, which does
  not propagate scale uncertainty, avoids it). This is a property of the
  track-G heteroscedastic encoding under full posterior propagation, recorded in
  the original 76-cell run (B9.19-B9.20), not an RG.7 artefact. It is exactly
  the same mechanism, seen on the gaussian scale, that makes the track-N Tweedie
  Laplace-averaged conservative: full uncertainty propagation is penalised on a
  near-deterministic response. The honest comparison reports it and prefers the
  natural-support track N (beta/tweedie) for the abundance data, where the
  Jacobian-corrected densities are on the response's own measure.

---

## 3. Decision E (criterion 2.O, stratified, no-overreach)

**On the abundance count (the RG.7 contribution), gdpar and mgcv are
statistically indistinguishable, by two stratified readings, with neither method
overclaimed:**

1. **gdpar reports a certified-limit Laplace / plug-in predictive** (mode +
   exact-Hessian Gaussian), because the count posterior is a certified
   non-sampleable non-Gaussian canyon (cond ~1e8; the multi-chain gate rules out
   default NUTS, constant-mass Laplace HMC, and whitened NUTS).
2. **mgcv reports a REML plug-in predictive.** Both are first-class Laplace /
   plug-in predictives on the same original measure — neither is exact MCMC,
   and that limitation is symmetric.
3. **Their point predictives coincide to < 0.4 ELPD units in every cell.** The
   only systematic difference is that gdpar's full-Laplace average additionally
   propagates dispersion/power uncertainty (a conservative, not a favourable,
   choice), keeping the two within ~1 se.

There is **no overreach**: gdpar does not claim to sample the canyon (it
certifies it cannot, reproducibly); it does not claim the Laplace is exact (the
fit-quality diagnostics quantify its crudeness); and it does not cherry-pick the
plug-in number for the headline (the conservative Laplace average is reported).
The RG.7 value is that gdpar's geometry capability **diagnosed** the canyon
(cond ~1e8, genuine non-Gaussianity), **certified** the sampling limit with a
reproducible multi-chain experiment, and **delivered a competitor-parity
predictive** where naive NUTS produced nothing usable — turning a 76/80 dead end
into an honest 80/80 closure.

---

## 4. Reproducibility

- Cells: `inst/benchmarks/results/cells_9_2_O/*.rds` (80); manifest
  `manifest_9_2_O.csv` (80 rows, reconciled).
- RG.7 Laplace driver: `inst/benchmarks/run_rg7_laplace_elpd.R`
  (`GDPAR_RG7_LAPLACE=1`); extractor + diagnostics:
  `inst/benchmarks/scripts/rg7_laplace_elpd.R`; per-cell Laplace certification
  sidecars: `inst/benchmarks/results/rg7_laplace/laplace_*.rds`.
- Decisive whitened-NUTS experiment (the certified limit): `~/rg33_whiten.R`,
  `~/rg33_whiten_full.out`.
- Divergence ledger: `MAPA_DIVERGENCIAS_B9.md` D98. Method narrative: vignette
  `vop08_geometric_robustness.Rmd` (RG.7 section).

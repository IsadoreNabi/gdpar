# Block 9.1 — Internal synthetic adversarial re-validation

**Status:** open (updated Session B9.8, 2026-05-28). Layers 1–3 substantive
+ drift smoke executed (B9.7). B9.8 was reoriented by user decision to build
**Axis 2 (dependence-robust inference, §10)** and resolve D74; the Stan-bound
layers of 9.1 (fit fuzz, 12-family outcome-outlier stress, full bit-exact
compare-path) remain gated and now move to the next session.

This report fulfils criterion 1 of the canonical closure of Block 9
(Charter §6). It is updated incrementally across the 9.1 sub-sessions.

---

## 1. Decision C — scope of the Path C fuzz protocol

Decision C (Charter §4), consulted batched at the start of B9.7 under
`feedback_upfront_session_consultations` with a six-dimensional analysis,
was resolved to **C.iv (lateral): layered defence with drift-gated
reseed** (canonized as **D73** in `MAPA_DIVERGENCIAS_B9.md`). C.iv
subsumes the three pure options:

- **Layer 1 (core of C.iii)** — structural adversarial fuzz over the
  canonical architecture: aliased grouping (C7), per-slot D-ID violations,
  degenerate priors, latent HOM stratification. Expected outcome: canonical
  aborts / `passed = FALSE` verdicts.
- **Layer 2 (core of C.i)** — property-based invariances: permutation and
  scale equivariance (deterministic, exact) plus noise-injection continuity
  and likelihood-exchangeability of the fitted estimates (statistical,
  Stan-bound, gated).
- **Layer 3 (bit-exact regression as a drift detector, subset of C.ii)** —
  re-run the frozen goldens against their RDS to confirm the cluster-K4
  Stan unification (9.3.a + cascade L + M.iv) preserved bit-exactness *by
  construction*.
- **Layer 4 (full reseed C.ii) — conditional** — triggered only if the
  drift detector in Layer 3 finds drift attributable to the unification.

The Charter-mandated invariants (extreme outcome outliers over the 12
canonical families; NA / missing robustness) are folded into Layers 1–3.

---

## 2. Pacing (Decision C.iv atomic sub-units, B9.7–B9.9)

Under `feedback_token_economy_session_planning`, the cheap-to-moderate
layers that need no heavy Stan compile/sampling were executed substantively
in B9.7; the Stan-bound layers are a separate nocturnal sub-unit (B9.8).

| Sub-unit | Content | Cost | Session |
|---|---|---|---|
| B9.7 | Layers 1–2 (deterministic) + 3 (NA/outlier algebra) test files + drift smoke (gaussian KxP) + report | fast (one ~1 s KxP fit) | **done** |
| B9.8 | Gated fit fuzz (INV-C) + 12-family outcome-outlier stress (STR-C) + full bit-exact compare-path Tier1 (~9 min) + KxP (~2.3 min) + EB goldens + Tier2 Tweedie (~38 min); conditional reseed if drift | Stan-heavy (nocturnal) | pending |
| B9.9 | Buffer + 9.1 closure + report finalization | — | pending |

---

## 3. Layer 1 — structural adversarial fuzz (`test-block9_adversarial.R`)

All aborts are pre-fit (input validation / identifiability algebra / prior
structure validation); the file runs by default with no Stan compilation.

**Naming correction (verified against the code):** the Charter names a
`gdpar_c7_violation` class; the actual canonical abort surface raises
`gdpar_input_error` or `gdpar_identifiability_error`, and the C4 / C4-bis
identifiability machinery returns a report with `passed = FALSE` rather
than aborting. The tests assert the *actual* behaviour.

| Test | Adversarial input | Asserted canonical outcome | Result |
|---|---|---|---|
| ADV-A1 | perfectly aliased additive covariate | `passed = FALSE`, collinear directions, κ > 1e6 | PASS |
| ADV-A2 | scaled-replica covariate (alias up to constant) | `passed = FALSE` (normalization cannot rescue) | PASS |
| ADV-A3 | degree-1 W aliased with additive block | per-coord C4-bis `passed = FALSE` | PASS |
| ADV-B1 | duplicate `prior_canonical_kind` across K≥2 slots | D-ID symbolic separability fails → report fails | PASS |
| ADV-B2 | clean K≥2 family (distinct kinds) | D-ID passes (negative control) | PASS |
| ADV-C1 | unknown family name | `match.arg` reject | PASS |
| ADV-C2 | the 12 canonical families | stan_id contract + canonical-kind tags hold | PASS |
| ADV-D1 | structurally degenerate prior args | `gdpar_input_error` | PASS |
| ADV-D2 | unknown / non-scalar `priors_by_kind` | `gdpar_input_error` | PASS |
| ADV-D3 | content-degenerate prior syntax (`normal(0,-1)`) | accepted R-side (documented boundary) | PASS |
| ADV-E1 | latent two-group HOM-violating slopes | pre-fit checks PASS (HOM is fit-time) | PASS |

**Boundary documented (ADV-D3):** `gdpar_prior()` validates prior
*structure* (non-empty character scalar, known canonical kind) but not the
numerical *content* of the Stan syntax; a negative scale is therefore
accepted at the R level and rejected only by the Stan parser/sampler. This
is the deliberate R-side / Stan-side validation boundary.

**Boundary documented (ADV-E1):** the homogeneity (HOM) assumption is a
distributional / fit-time regularity condition, not a pre-fit algebraic
invariant. Latent-group heterogeneity leaves the design full rank, so the
pre-fit identifiability check passes; the violation can only surface
post-fit (e.g. via `gdpar_contraction_diagnostic`). A post-fit HOM stress
test is a candidate for the B9.8 gated layer.

---

## 4. Layer 2 — property-based invariances (`test-block9_invariances.R`)

| Test | Invariance | Level | Result |
|---|---|---|---|
| INV-A1 | row permutation ⇒ identifiability eigenstructure invariant | deterministic | PASS |
| INV-A2 | row permutation ⇒ per-coord C4-bis verdict invariant (p>1) | deterministic | PASS |
| INV-B1 | covariate rescaling ⇒ verdict + conditioning invariant | deterministic | PASS |
| INV-B2 | covariate sign flip ⇒ eigenstructure invariant | deterministic | PASS |
| INV-C1 | EB location estimate permutation-invariant (MC error) | statistical (GATED) | deferred B9.8 |
| INV-C2 | EB location estimate continuous under noise injection | statistical (GATED) | deferred B9.8 |

The deterministic invariances hold exactly (up to floating-point summation
order): the Gram matrix `Z'Z` is a sum over rows (permutation-invariant)
and the check normalizes columns to unit norm (scale/sign-invariant).

---

## 5. Layer 3 — stress outliers + NA robustness (`test-block9_stress.R`)

| Test | Stress | Asserted outcome | Result |
|---|---|---|---|
| STR-A1 | NA in an AMM-spec covariate | `gdpar_input_error` (pre-fit) | PASS |
| STR-A2 | NA in the outcome | `gdpar_input_error` (clean abort) | PASS |
| STR-A3 | NaN in the outcome | `gdpar_input_error` (treated as missing) | PASS |
| STR-B1 | 99.9% / 0.1% quantile covariate outliers | finite conditioning, full rank | PASS |
| STR-B2 | covariate on a 1e6 scale | finite conditioning (normalization) | PASS |
| STR-C  | outcome outliers over all 12 families | finite estimates ∨ canonical abort | GATED B9.8 |

---

## 6. Layer 3 (drift detector) + Layer 4 (conditional reseed)

A fast bit-exact drift smoke was run in B9.7 on the freshest golden of the
cluster-K4 closure: **`golden_fb_KxP_gaussian_polynomial_K2_p2`** (seed
91001, produced in B9.6). The gaussian K=2, p=2 FB fit was reproduced in
isolation (without overwriting the golden) and compared against the frozen
RDS:

```
DRIFT>> bit-exact identical: TRUE
DRIFT>> fresh dim: 400x815 | golden dim: 400x815
```

**Verdict:** no drift. The cluster-K4 Stan unification (9.3.a dispatcher +
cascade L EB-side + M.iv unified assembler + N.a fit harness) preserved
bit-exactness *by construction*, as asserted across B9.3–B9.6. Therefore
**Layer 4 (full reseed C.ii) is NOT triggered**; a full reseed would be
redundant compute.

The full bit-exact compare-path (Tier1 12 configs ~9 min + KxP 3 remaining
seeds + Tier2 Tweedie K=3 ~38 min + EB 8.6.B/C/D goldens) is the Stan-heavy
B9.8 sub-unit and remains the comprehensive confirmation; the B9.7 smoke is
the leading indicator.

---

## 7. Findings — candidate robustness refinements (D74; finding 1 resolved in B9.8, see §10.2)

The adversarial protocol surfaced two behaviours that meet the Charter's
"abort clean, no silently-incorrect results" bar but expose a refinement
opportunity. They are **documented, not silently patched** (the fix touches
~6 shared validation sites; reserved for an upfront B9.8 consultation under
`feedback_nocturnal_autonomy`).

1. **Non-finite (Inf/-Inf) outcome is not pre-validated.** The outcome
   check uses `is.na()` only; `Inf` in `y` passes the pre-fit guard and
   surfaces downstream as a `gdpar_eb_numerical_error` (Stan failure)
   rather than a clean `gdpar_input_error`. Candidate root-cause fix:
   extend the `is.na(y)` guard to `is.na(y) | !is.finite(y)`, symmetric to
   the existing NA check, at the ~6 outcome-validation sites
   (`R/gdpar.R`, `R/eb.R`). `NaN` is already caught (it is `is.na`).
2. **NA in a formula-RHS-only covariate is dropped, not aborted.** NA in
   an AMM-spec covariate (a/b/W/x_vars) aborts (`R/amm_spec.R`), but NA in
   a covariate that appears only in the model `formula` RHS is dropped by
   R's `model.frame` `na.action` default. This is standard, documented R
   behaviour (not silently-incorrect), but it is asymmetric with the
   AMM-covariate guard.

Both are **candidate** refinements, not defects: in every case the pipeline
either aborts or follows documented R semantics; none silently returns an
incorrect fit.

---

## 8. New test inventory (B9.7)

- `tests/testthat/test-block9_adversarial.R` — 11 test_that (Layer 1).
- `tests/testthat/test-block9_invariances.R` — 6 test_that (Layer 2; 4
  default + 2 gated `GDPAR_RUN_BLOCK9_FIT_FUZZ`).
- `tests/testthat/test-block9_stress.R` — 6 test_that (Layer 3; 5 default +
  1 gated `GDPAR_RUN_BLOCK9_STRESS_FITS`, 12-family roster).

All default-running tests PASS; gated tests SKIP by default.

---

## 9. Criterion-1 checklist (Charter §6.1)

| Component | B9.7 | B9.8 | B9.9 (Stan-bound) |
|---|---|---|---|
| Fuzz protocol Path C (adversarial structural) | ✅ executed | — | — |
| Identifiability adversarial (C7 / D-ID / HOM) | ✅ executed | — | — |
| Property-based invariances | ✅ deterministic | — | ✅ statistical fit fuzz (INV-C1/C2) |
| Bit-exact re-bootstrap (drift) | ✅ smoke (KxP gaussian) | — | ✅ full Tier1/Tier2/KxP/EB (no drift) |
| Stress outliers (99.9% / 0.1%) | ✅ covariate-side | — | ✅ outcome-side 12 families (STR-C) |
| NA / missing robustness | ✅ executed | Inf guard (D74) | — |
| Reseed (conditional C.ii) | not triggered (no drift) | — | not triggered (no drift; §11) |

Criterion 1 is **SATISFIED and CLOSED (B9.9).** B9.7 executed the
deterministic / cheap layers; B9.8 added Axis 2 (dependence-robust inference,
§10) ahead of the Stan-bound work; B9.9 ran the deferred Stan-bound layers
(§11) with zero drift and zero test failures.

---

## 10. B9.8 — dependence-robust inference (Axis 2) + D74 resolution

### 10.1 Scope reorientation (user decision)

At the start of B9.8 the user reopened the autocorrelation question. The
6-dimensional analysis separated two axes:

- **Axis 1 — modelling the dependence** (autoregressive errors, spatial
  random fields). Architecturally a new modelling block of Block-8 size;
  computationally graded (AR(1) is cheap in Stan, continuous-domain SPDE is
  INLA's home turf); deferred to a future block (Block 10), evidence-gated
  by the external validation 9.2. The honest weakness stands: gdpar does
  not model residual correlated noise or spatial random effects.
- **Axis 2 — valid inference *under* dependence** without modelling it.
  In scope now under the robustness-until-B9 regime; built this session.

A code audit established the decisive fact: gdpar's only pre-existing "block
bootstrap" (`R/preflight.R`) resamples *MCMC draws within a chain* for the
CP/NCP info-ratio test — it is an MCMC-autocorrelation device, **not** a
data-dependence-robust SE estimator. The reported uncertainty (the
posterior, EB and FB alike) is therefore not robust to data dependence on
any path. Axis 2 is net-new machinery, not a tune-up. Consequently the
Stan-bound layers of §9 were moved to the next session and B9.8 delivered
the Axis 2 atom below.

### 10.2 D74 resolution (finding 1)

Finding 1 of §7 is **resolved**: the six outcome-validation sites in
`R/gdpar.R` (3) and `R/eb.R` (3) now use a type-aware non-finite predicate
(`!is.finite(y)` for numeric outcomes, `is.na(y)` otherwise). `Inf` / `-Inf`
outcomes abort pre-fit with `gdpar_input_error` (previously a late, opaque
`gdpar_eb_numerical_error`). `test-block9_stress.R` gains STR-A4. Finding 2
(formula-RHS-only NA) remains deferred (a design choice, not a one-liner).

### 10.3 The diagnostic (D75)

`gdpar_dependence_diagnostic()` measures serial dependence in the residuals
of a scalar EB fit (lag-1 autocorrelation with a normal-approximation
p-value, Durbin-Watson, Ljung-Box). It converts the iid-violation hazard
into a measured quantity and gates the remedy. Smoke validation on injected
AR(1) noise (rho = 0.7, n = 120): lag-1 autocorrelation 0.68 (p ~ 9e-14),
Durbin-Watson 0.595, Ljung-Box Q(20) = 159 (p ~ 0) -> "Dependence
detected"; the matched iid control did not flag.

### 10.4 The remedy (D75)

`gdpar_dependence_robust()` reports dependence-robust SEs and percentile
intervals via a temporal moving / circular block bootstrap that refits the
model on contiguous-block resamples (default block length n^(1/3); Künsch
1989). This is the working-independence + robust-variance stance of Liang &
Zeger (1986), separated by rigor layer:

- **Algebraic.** The estimand is the EB coefficient vector
  {theta_ref, a_coef, b_coef, W_raw}; the resampler draws contiguous blocks
  ordered by the user-supplied temporal index.
- **Statistical.** Valid for weak / short-range dependence relative to the
  block length; it does not rescue long-memory or unit-root processes.
  Point estimates are unchanged (consistent under a correct mean structure,
  not efficient); only the uncertainty is made robust.
- **Numerical.** Per-refit Stan seeds are derived deterministically from the
  master seed; same-seed runs reproduce the robust SEs on a fixed toolchain.

### 10.5 Golden decision (D76)

The bit-exact golden is the **deterministic resampler** (pure R RNG, frozen
inline in `test-block9_dependence.R`). The HMC-dependent bootstrap SEs are
**not** frozen bit-exact: HMC is not reproducible across Stan toolchains, so
a frozen fit-based golden would be a fragile, masking fixture contrary to
the no-check-masking regime. The HMC-dependent behaviour is validated by the
gated end-to-end tests instead.

### 10.6 Test inventory (B9.8)

- `tests/testthat/test-block9_dependence.R` — Stan-free resampler algebra +
  API guards + resampler golden (default-running); AR(1) detection +
  robust-SE table structure + seed-determinism gated by
  `GDPAR_RUN_BLOCK9_DEP_FITS`.
- `test-block9_stress.R` STR-A4 (Inf outcome aborts pre-fit; D74).

### 10.7 Deferred to the next session

Spatial diagnostic (Moran's I) + spatial block bootstrap; the full-Bayes
path; Politis-White (2004) automatic block length; K > 1 / p > 1; and the
Stan-bound layers of §9 (fit fuzz, 12-family outcome-outlier stress, full
compare-path drift detection).

### 10.8 Validation results (B9.8 close)

- **Default suite** (Stan-free): PASS (17 PASS / 2 SKIP gated in
  `test-block9_dependence.R`; resampler algebra + API guards + resampler
  golden D76).
- **Gated Layer 3** (`GDPAR_RUN_BLOCK9_DEP_FITS=1`): PASS. The diagnostic
  flags AR(1) and not iid; `gdpar_dependence_robust()` returns a finite
  robust-SE table with ordered percentile intervals and is seed-deterministic
  across same-seed runs. (8 benign HMC divergences at 0.0% from the tiny
  example-grade iterations.)
- **R CMD check natural `--as-cran --run-donttest`** (no
  `_R_CHECK_FORCE_SUGGESTS_` masking): **Status 2 NOTEs**, 0 ERROR, 0 WARNING
  — identical to the S22 / B9.5 / B9.6 / B9.7 baseline (CRAN incoming
  feasibility dev version + examples timing D59). Tests [15m] OK; Rd and
  code/documentation checks OK (no new NOTE from the B9.8 code). A first
  run surfaced one cosmetic test regression (`test-gdpar_multi.R` asserted the
  old "missing" wording of the outcome guard, now "non-finite"); the regexp
  was updated and the `\donttest` examples were lightened (B = 10), after
  which the natural check returned to the clean baseline. Logs:
  `~/claude-setup/gdpar-block-9/logs/B9_8_gated.log` and
  `B9_8_rcmdcheck_post.log`.

---

## 11. B9.9 — Stan-bound layer closure (Criterion 1 closed)

The layers B9.8 deferred were run as one nocturnal sub-unit
(`~/claude-setup/gdpar-block-9/logs/run_B9_9_gated.R`, log
`B9_9_gated.log`, ~59 min wall clock). The drift detector is a bit-exact
re-fit compared against the frozen RDS; any mismatch or test failure is the
Layer-4 escalation signal (reseed C.ii + user consult). The canonical
goldens were never overwritten: the KxP fit closures return the fit object,
and the EB re-fits were redirected to a temp directory. Order put the cheap,
most-diagnostic compares first and the ~38 min Tweedie K=3 Tier2 last.

### 11.1 Bit-exact compare-path (drift detector) — no drift

| Path | Configs | Result |
|---|---|---|
| KxP FB Path C (re-fit closures) | g_KxP_002 beta K2 p2 (91002); g_KxP_003 gamma K2 p3 (91003); g_KxP_004 student_t K3 p2 (91004) | **bit-exact identical = TRUE** (400×815 / 400×1521 / 400×980); g_KxP_001 gaussian was the B9.7 smoke |
| EB 8.6.B scalar (re-fit into temp dir) | gaussian/poly, poisson/poly, neg_binomial_2/poly, bernoulli/poly, gaussian/bspline | **draws bit-exact = TRUE and `theta_ref_hat` exact = TRUE** for all 5 |
| 8.3.9 Tier1 (`GDPAR_GOLDEN_8_3_9_COMPARE_TIER1=1`) | 12 configs (~9 min) | **PASS** (bit-exact) |
| 8.3.9 Tier2 (`GDPAR_GOLDEN_8_3_9_COMPARE_TIER2=1`) | Tweedie K=3 (~38 min) + bspline_p2 | **PASS** (bit-exact) |

**Honest scope note (EB).** Only the 5 fitable 8.6.B EB goldens carry a
`draws_matrix` on disk and are therefore compare-able. The 8.6.C (13
fitable) and 8.6.D (4 fitable) entries were never persisted as fitable RDS
(only their guard / metadata entries exist; the fitable bootstrap was
deferred, consistent with the D58 tarball-size posture), so there is no
frozen artefact to drift against — they are out of scope for a bit-exact
compare, not skipped silently.

**Verdict:** no drift on any of the 22 compare-able goldens (3 KxP + 5 EB +
12 Tier1 + 2 Tier2). The cluster-K4 Stan unification (9.3.a dispatcher +
cascade L EB-side + M.iv unified assembler + N.a fit harness) preserved
bit-exactness, as asserted by construction across B9.3–B9.6 and confirmed
comprehensively here. **Layer 4 (full reseed C.ii) is NOT triggered.**

### 11.2 Gated package tests — `[ FAIL 0 | WARN 8 | SKIP 0 | PASS 76 ]`

- **INV-C1 / INV-C2** (`GDPAR_RUN_BLOCK9_FIT_FUZZ=1`): PASS. The EB location
  estimate is permutation-invariant up to MC error and Lipschitz-continuous
  under O(1e-3) noise injection.
- **STR-C** (`GDPAR_RUN_BLOCK9_STRESS_FITS=1`): PASS. All 12 canonical
  families fed outcome outliers at the 99.9 % / 0.1 % tails returned finite
  `theta_ref_hat` or aborted with a canonical error; none silently produced
  NaN/Inf estimates.
- `SKIP = 0` proves every gated `test_that` (including both compare-path
  tiers) executed its assertions — nothing skipped on a missing RDS or unset
  gate.

### 11.3 D77 — candidate refinement surfaced by the first substantive STR-C run

The 8 warnings are all from the gated STR-C fits (`test-block9_stress.R:138`)
and all of one kind: `compute_diagnostics()` (`R/gdpar.R:1724`) computes
`max(rhat_vec, na.rm = TRUE)`; under degenerate outlier-laden data the whole
R-hat vector can come back all-`NA`, so `max()` warns
("no non-missing arguments to max; returning -Inf") and the reported
`max_rhat` diagnostic becomes `-Inf`. This is the first substantive run of
STR-C (deferred B9.7 → B9.8 → B9.9), so it is a genuinely new observation.

It meets the Charter bar (no silently-incorrect result: the point estimate
stays finite — the test PASSES — and `-Inf` plus a warning is a visible flag,
not a quiet corruption) but exposes a refinement opportunity, exactly like
the D74 findings. Candidate root-cause fix: guard the all-`NA` case
(`if (all(is.na(rhat_vec))) NA_real_ else max(rhat_vec, na.rm = TRUE)`).
**Documented, not patched** (the patch touches shared diagnostic code; under
nocturnal autonomy and the D74 precedent it is reserved for an upfront
consultation). Confined to the gated stress path — it never fires in the
default suite or in `R CMD check`, so the baseline is unaffected. Logged as
**D77** in `MAPA_DIVERGENCIAS_B9.md` §2.1.

### 11.4 Criterion 1 — CLOSED

Every component of the Charter §6.1 / §2.2 criterion-1 checklist is now
green (§9). Internal synthetic adversarial re-validation is complete: the
adversarial / identifiability / invariance / stress / NA layers (B9.7), the
dependence-robust inference axis (B9.8, §10), and the Stan-bound fit-fuzz /
12-family stress / full bit-exact compare-path (B9.9) all pass with no drift
and no failures. The one candidate refinement (D77) is documented and
deferred. **Sub-bloque 9.1 is closed; Bloque 9 proceeds to 9.2 (external
re-validation vs TOP3).**

### 11.5 Natural `R CMD check` (no masking)

`R CMD check --as-cran --run-donttest` on the freshly built tarball, with no
`_R_CHECK_FORCE_SUGGESTS_` masking (`feedback-no-check-masking-b9`):
**Status 2 NOTEs, 0 ERROR, 0 WARNING (exit 0)** — identical to the
S22 / B9.5 / B9.6 / B9.7 / B9.8 baseline. NOTE 1 is the CRAN incoming
feasibility note (new submission + development version `0.0.0.9001` +
`cmdstanr` / `INLA` suggests resolved via `Additional_repositories`); NOTE 2
is the examples-timing note (D59). Examples [26m/27m], tests, and vignette
rebuild [205s] all OK. B9.9 introduced no new NOTE/WARNING/ERROR — expected,
since it is a validation-only session with no source changes. Log:
`~/claude-setup/gdpar-block-9/logs/B9_9_rcmdcheck_post.log`.

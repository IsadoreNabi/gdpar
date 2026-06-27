# Block RG / RG.1.c -- posterior-geometry diagnostic calibration

Generated 2026-06-02 by `inst/benchmarks/rg_calibration_analysis.R` over the cells in `cells_rg`.

This report replaces the **8/8 diagonal of one replica** smoke of B9.22 with **measured error rates**: the `gdpar_geometry_diagnostic` classifier run over a `gdpar_geometry_suite` grid of three axes -- difficulty (severity ladder), pilot budget (cheap 100/100 vs default 150/150), and replicas -- with an adaptive, minimax-robust replica allocation and an out-of-sample (held-out) threshold calibration.

## Methodology

- **Signals are collected once** by cheap NUTS pilots through cmdstan; the threshold sweep is a pure function of the stored signals (no resampling). Size-invariant signals only -- divergences, E-BFMI, tree depth, condition number, step/scale -- **never R-hat/ESS on short runs** (lesson B9.20/B9.21).
- **Adaptive allocation (Q2 mix 3+4):** Stage 1 gives every cell_base R0=5 replicas; Stage 2 tops up to R=25 only the cell_bases whose Wilson confidence interval for the correct rate is *unresolved* (overlaps 0.5, or rate 0). Minimax-robust: replicas are spent where the verdict is statistically undecided; settled cells keep R0.
- **Threshold calibration (Q3 union):** a data-driven Youden-optimal cut per gating signal *proposes* thresholds; a coordinate-ascent grid around the interpretable defaults *regularises* them on a CALIBRATION fold (odd replicas); the final thresholds are *reported* on a held-out VALIDATION fold (even replicas). The classifier stays the transparent rule-based one; only its numeric thresholds move.
- **Honesty (no overreach):** these are measured rates, not a curated diagonal. gdpar's diagnostic *can* and *does* fail targets -- mild pathologies look benign, and some pathologies are mutually confusable. Reported below as-is.

## Coverage: 540 cells (240 Stage 1 + 300 Stage 2), 48 cell_bases.

## Headline accuracy

**Design-uniform (macro) accuracy** is the honest population estimate: each cell_base (target x difficulty x pilot) is weighted equally, so the adaptive Stage-2 top-up (which gives hard cells more replicas) does not bias it. The replica-pooled figures are the tuning objective and are reported alongside but are pessimistic on the adaptive sample (they over-weight the hard, over-sampled cells).

| thresholds | macro accuracy (design-uniform) | pooled accuracy | pooled balanced |
|---|---|---|---|
| default | 0.781 | 0.513 | 0.625 |
| tuned | 0.908 | 0.826 | 0.901 |

Out-of-sample (held-out validation fold) balanced accuracy -- the overfitting guard:

| thresholds | held-out accuracy | held-out balanced |
|---|---|---|
| default | 0.463 | 0.598 |
| tuned | 0.805 | 0.894 |

Calibration-fold balanced accuracy at the tuned thresholds: 0.912.

Design-uniform per-class accuracy (tuned):

`anisotropic=0.93  boundary=1.00  flat_direction=1.00  funnel=0.63  heavy_tails=0.71  isotropic=1.00  multimodal=1.00  quasi_deterministic=1.00`

## Held-out confusion matrix (tuned thresholds)

| truth \\ pred | isotr | aniso | funne | heavy | quasi | multi | bound | flat_ |
|---|---|---|---|---|---|---|---|---|
| isotropic | 12 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| anisotropic | 0 | 32 | 0 | 0 | 0 | 0 | 0 | 0 |
| funnel | 0 | 0 | 39 | 2 | 5 | 11 | 0 | 5 |
| heavy_tails | 0 | 2 | 17 | 27 | 6 | 0 | 0 | 0 |
| quasi_deterministic | 0 | 0 | 0 | 0 | 12 | 0 | 0 | 0 |
| multimodal | 0 | 0 | 0 | 0 | 0 | 12 | 0 | 0 |
| boundary | 0 | 0 | 0 | 0 | 0 | 0 | 32 | 0 |
| flat_direction | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 32 |

Class abbreviations are the first five letters (isotr, aniso, funne, heavy, quasi, multi, bound, flat_).

## Held-out per-class recall (tuned thresholds), Wilson 95% CI

| class | n | correct | recall | ci_lo | ci_hi |
|---|---|---|---|---|---|
| isotropic | 12 | 12 | 1.000 | 0.758 | 1.000 |
| anisotropic | 32 | 32 | 1.000 | 0.893 | 1.000 |
| funnel | 62 | 39 | 0.629 | 0.505 | 0.738 |
| heavy_tails | 52 | 27 | 0.519 | 0.387 | 0.649 |
| quasi_deterministic | 12 | 12 | 1.000 | 0.758 | 1.000 |
| multimodal | 12 | 12 | 1.000 | 0.758 | 1.000 |
| boundary | 32 | 32 | 1.000 | 0.893 | 1.000 |
| flat_direction | 32 | 32 | 1.000 | 0.893 | 1.000 |

## Error-vs-severity (accuracy by ground-truth class x difficulty)

Difficulty rank low/mid/high = mild/default/severe pathology intensity (mid ~ the suite default).

| class | low | mid | high |
|---|---|---|---|
| anisotropic | 1 | 1 | 0.8 |
| boundary | 1 | 1 | 1 |
| flat_direction | 1 | 1 | 1 |
| funnel | 0.84 | 0.97 | 0.14 |
| heavy_tails | 1 | 0.7 | 0.38 |
| isotropic | 1 | 1 | 1 |
| multimodal | 1 | 1 | 1 |
| quasi_deterministic | 1 | 1 | 1 |

## Accuracy-vs-pilot-cost (tuned thresholds)

| pilot budget | accuracy | balanced | n |
|---|---|---|---|
| cheap | 0.814 | 0.896 | 280 |
| default | 0.838 | 0.900 | 260 |

## Thresholds: default vs tuned

| threshold | default | tuned | changed |
|---|---|---|---|
| divergent_rate_high | 0.01 | 0.01 |  |
| funnel_ebfmi_low | 0.25 | 0.35 | yes |
| heavy_cond_max | 8 | 25 | yes |
| treedepth_sat_high | 0.2 | 0.2 |  |
| condition_high | 50 | 12 | yes |
| step_scale_ratio_low | 0.1 | 0.1 |  |
| nslope_grows | 0.3 | 0.8 | yes |
| flat_var_high | 1000 | 600 | yes |
| boundary_prox_high | 0.1 | 0.02 | yes |
| boundary_eps | 0.01 | 0.01 |  |
| multimodal_high | 2 | 2.5 | yes |
| heavy_kurtosis_high | 3 | 1.8 | yes |
| target_ess | 400 | 400 |  |

## Anchors (no overreach)

Girolami & Calderhead (2011, RMHMC); Betancourt (SoftAbs, E-BFMI, pathology taxonomy); Neal (the funnel); Lu et al. (relativistic Monte Carlo); Pennec / Bhattacharya-Patrangenaru (Riemannian statistics, the user's `ON_RIEMMANIAN_STATISTICS`). gdpar integrates and orchestrates these with its structural knowledge; it does not claim to invent the methods, nor that its classifier is unique or infallible.


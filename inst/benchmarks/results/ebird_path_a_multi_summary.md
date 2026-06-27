# Block 7 - eBird path-A multi: verdict summary

Generated at: 2026-05-19 19:17:00

## Scenario

- Single scenario `ebird_path_a_multi` (p = 2: logit-encounter / log1p-count).
- 5 species (reevir1, amered, ovenbi1, btnwar, woothr) x
  4 NE-USA breeding sub-regions = 20 groups.
- 12000 train / 3000 test rows after stratified sub-sampling.
- Predictors in `a`: lon, lat, bio1, bio12, bio15 (WGS84).
- Ground truth: S&T seasonal-breeding rasters
  (y1 = logit(occurrence_mean), y2 = log1p(count_mean)).

## Ranking by elpd_loo (higher is better)

| rank | method | elpd_loo | se | delta vs best |
|------|--------|----------|----|---------------|
| 1 | mgcv | -28997.12 | 188.23 | 0.00 |
| 2 | rstanarm | -28999.24 | 188.28 | -2.12 |
| 3 | gdpar | -28999.74 | 188.36 | -2.62 |
| 4 | brms | -28999.74 | 188.28 | -2.63 |
| 5 | inla | -29000.34 | 188.41 | -3.22 |

## Wall time per method (seconds)

| method | wall_seconds |
|--------|--------------|
| gdpar | 6600.30 |
| brms | 1126.16 |
| mgcv | 0.71 |
| inla | 2.82 |
| rstanarm | 244.11 |

## Verdict

- Best elpd_loo: **mgcv** (-28997.12, se 188.23); runner-up: **rstanarm**.
- gdpar rank: **3 / 5**.
- elpd_loo spread across all 5 methods: 3.22 (reference se ~188.3); differences are within 1 se, i.e. statistically indistinguishable on this dataset.
- PSIS-LOO health: all methods Pareto-k < 0.7 (PSIS-LOO healthy).

## Caveats

- **Model-vs-model validation, not ground-truth validation.** Both fitted targets (y1, y2) are themselves model outputs from the Cornell Lab Status & Trends pipeline (an ensemble of boosted trees), not direct field observations. The bench compares which framework recovers a smoothed S&T surface best, not which one approximates actual occurrence/count.
- **Sub-region fallback is state-based, not BCR-based.** BCRs (Bird Conservation Regions) are the ecologically principled grouping for North American passerines; this run uses a state-based fallback (NE_FOREST = ME+NH+VT, MID_ATL_COAST = MA+RI+CT+NJ+DE+MD, APPALACHIAN = PA+WV, GREAT_LAKES = NY). NY entire is the worst offender: it crosses 4 BCRs and is grouped as a single sub-region, which dilutes the random effect.
- **No competitor was tuned.** All 5 frameworks ran with default hyperpriors / smoother bases / mesh defaults. Differences in absolute elpd_loo are tiny (within noise) so the ranking should not be over-interpreted.
- **gdpar wall is dominated by Stan compilation + 4 chains x 4 cores HMC**; mgcv / INLA are GLMM analogues with PIRLS / INLA approximation. Wall comparison is informative for engineering, not for statistical efficiency.


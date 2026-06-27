# gdpar

**General Dynamic Parameter Models via Reference Anchoring**

`gdpar` is an R package implementing a unified predictive framework in which each
unit's parameter is decomposed around a population reference:

```
theta_i = theta_ref + Delta(x_i, theta_ref)
```

with `Delta` following the canonical **Additive--Multiplicative--Modulated (AMM)**
form `a(x) + b(x) * theta + W(theta) x`. Here `theta_ref` is a population
reference estimated from the data and `Delta` is an explicit, identifiable
deviation function. The package provides identifiability diagnostics, validity
tests for the population reference, distributional (multi-slot) regression,
empirical-Bayes and full-Bayes estimation, a causal-inference bridge, a
geometric-robustness sampling engine, and dependence-robust inference.

## Installation

Path 1 (the operational estimation path) uses Stan via `cmdstanr`. Install the
back-end once per machine:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

Then install `gdpar` from GitHub:

```r
# install.packages("remotes")
remotes::install_github("IsadoreNabi/gdpar")
```

## Quick start

```r
library(gdpar)

set.seed(20260526)
n <- 200L
x <- rnorm(n)
y <- 1.5 + 0.8 * x + rnorm(n, sd = 0.3)   # theta_ref = 1.5, a = 0.8
dat <- data.frame(y = y, x = x)

fit <- gdpar(
  formula = y ~ x,
  data    = dat,
  family  = gdpar_family("gaussian"),
  path    = "bayes"
)

summary(fit)
coef(fit)
predict(fit, newdata = data.frame(x = seq(-2, 2, length.out = 11)), level = 0.9)
```

See the `vop00_quickstart` vignette for a line-by-line walk-through.

## The three estimation paths

`gdpar` specifies three estimation paths for the AMM decomposition:

- **Path 1 --- hierarchical Bayesian** (`path = "bayes"`): the operational path,
  via Stan / `cmdstanr`.
- **Path 2 --- varying-coefficient** (`path = "vcm"`): reference-grade asymptotic
  theory; implementation queued for a future version.
- **Path 3 --- hypernetwork** (`path = "hyper"`): reference-grade; implementation
  queued for a future version.

Calls to the not-yet-implemented paths abort with
`gdpar_unsupported_feature_error`.

## Capabilities

- **Distributional regression (`K >= 1` slots):** Gaussian, Student-t, Gamma,
  Beta, Tweedie, Poisson, negative-binomial, zero-inflated and hurdle families,
  plus heterogeneous per-slot families; B-spline `W` bases; arbitrary covariate
  dimension `p`; grouped references.
- **Empirical-Bayes and full-Bayes** estimation, with EB-vs-FB comparators
  (marginal total variation and joint KSD).
- **Causal positioning:** a T-learner CATE/ITE bridge on the AMM side and an
  external meta-learner comparator (grf, EconML via reticulate).
- **Geometric-robustness sampling engine:** an opt-in, default-bit-identical
  Riemannian / sub-Riemannian / Finsler metric stack with a certifying
  orchestrator, for stiff and near-deterministic posteriors.
- **Dependence-robust inference:** diagnostics (lag-1 / Durbin-Watson; Moran's I)
  and block-bootstrap robust standard errors for temporally and spatially
  dependent data, with a data-driven (Politis-White) block length, on both the
  empirical-Bayes and full-Bayes paths. `gdpar` does not model the dependence;
  only its inference is made robust to it.
- **Residual diagnostics** (Dunn-Smyth quantile residuals; G1/G2/G3) and a
  golden-regression test suite.

## Validation

`gdpar` was re-validated against `mgcv`, `brms`, `INLA` and `rstanarm`:

- A **synthetic recovery** benchmark with known ground truth (3200 cells):
  `gdpar` leads on the distributional structure it models (heavy-tail quantiles,
  zero-inflated count means, heteroscedastic scale), is robust where it makes no
  modelling claim (e.g. autocorrelation it does not model), and is the most
  computationally expensive method in the roster.
- An **organic eBird** benchmark (80 cells): `gdpar` is statistically
  indistinguishable from `mgcv` on the held-out predictive density.

Reports: `inst/benchmarks/results/block9_synthetic_recovery.md` and
`inst/benchmarks/results/block9_revalidation.md`.

## Documentation

The package ships an extensive set of vignettes: an operational track
(`vop00`--`vop09`, from quickstart to dependence-robust inference) and a
theoretical track (`v00`--`v09`, covering the canonical form, identifiability,
the gnoseological validity of the population reference, the per-path asymptotic
theory, the empirical-Bayes versus full-Bayes comparison, and the
causal-inference positioning).

## Author and citation

José Mauricio Gómez Julián (ORCID
[0009-0000-2412-3150](https://orcid.org/0009-0000-2412-3150)).

## License

GPL (>= 3).

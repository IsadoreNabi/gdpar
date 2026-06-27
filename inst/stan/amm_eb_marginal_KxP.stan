// gdpar Path 1 -- Empirical Bayes (EB) Marginal Template, Path C
// (Sub-phase 8.6.D, decisions D36 = alpha + D37 = (i) canonized in
// Session 13a / 2026-05-25; decision D38'' = (h) canonized in
// Session 13b / 2026-05-25 after lateral analysis revealed a deadlock
// between the (b') K-parallel-outcome interpretation and the existing
// gdpar_formula_set single-outcome API).
//
// Architecture (D38'' = (h) canonical Path C, Session 13b):
//
//   Composition of Path A multivariate (Theorem 7A* of v07b) with
//   Path B multi-parametric (Theorem 7C* of v07b): single outcome
//   y[n, p] multivariate (Path A coord-wise factorization across p),
//   K slots distributional regression per coord (Path B multi-
//   parametric canonical link per slot).
//
//   - Outcome y_real[n, p] / y_int[n, p]: matrix-column outcome
//     shared across the K slots, coord-wise factorized across p.
//     The K slots are the canonical distributional parameters of the
//     family (e.g., (mu, log sigma) for Gaussian K=2, (log mu, log
//     phi) for NB K=2). Each slot has its own anchor in R^p per
//     group, mirroring Path A.
//   - theta_ref_kp[J_groups, K] vector[p]: per-group, per-slot,
//     per-coord anchor. Total J_groups * K * p anchor values.
//   - Likelihood (coord-wise factorization, multi-parametric family
//     per coord):
//        L = prod_{i=1}^n prod_{j=1}^p
//            D_{family_id_k[1]}(y[i, j] | eta_kp[1, i, j],
//                                          inv_link_2(eta_kp[2, i, j]),
//                                          ...)
//     where eta_kp[k, i, j] is the linear predictor of slot k at
//     observation i and coord j.
//
// Coverage (decision D40' = direct execution of D36 narrowed by D38''
// canonization to multi-parametric Path B-set):
//
//   Initial iteration of 8.6.D (Session 13b): stan_id in {1, 3}
//   (Gaussian K=2 and Negative Binomial K=2). These are the two
//   K=2 families with the most robust Path B numerical profile and
//   match the "Gaussian/NB confirmadas" reading of D36 when
//   interpreted as their K=2 versions in the Path B set
//   {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13}.
//
//   Remaining Path B families {5 Beta, 6 Gamma, 7 Lognormal_loc_scale,
//   8 Student-t, 9 Tweedie, 10 ZIP, 11 ZINB, 12 Hurdle-Poisson,
//   13 Hurdle-NB} deferred to a later iteration of 8.6.D with the
//   explicit numerical caveat of opening Section 6.1 (HMC condicional
//   bajo plug-in EB cerca del borde de soporte logit/log links + warmup
//   corto). Heterogeneous slot families (sub-phase 8.3.7 pattern)
//   deferred to a later iteration.
//
// Restrictions (registered as explicit debts against Block 9.x):
//
//   - use_W = 0 enforced (D39): the modulating component W is
//     orthogonal to the K x p cross structure and its canonical
//     factorization across both axes is registered as a follow-on
//     debt. Path C first iteration is pure (a, b, theta_ref) plus
//     per-slot K = 2 dispersion / power.
//   - a / b coefficients are per-slot per-coord ragged (D41): each
//     (slot, coord) pair has its own J_a_per_kp[k, j] /
//     J_b_per_kp[k, j] free columns. Maximum statistical generality.
//
// Family identifiers per slot (K=2 distributional regression families
// of the Path B set, applied coord-wise across p coords of the
// outcome):
//   1 = gaussian K=2          (mu identity, sigma log) -> per coord
//   3 = neg_binomial_2 K=2    (mu log,      phi log)   -> per coord
//
// Internal sampling parametrization (preserved from Path B K):
//   a : NCP only per slot (per-slot scale sigma_a_k; CP / per-coord
//       per-slot deferred).
//   b : LINEAR reparametrization per slot per coord on c_b_kp =
//       theta_ref_kp * b_coef_kp; NCP sample of c_b_kp_raw and
//       derive b_coef_kp post-hoc in the single-anchor regime.
//   theta_ref_kp : per-group, per-slot, per-coord anchor.
//
// The companion conditional template amm_eb_conditional_KxP.stan
// moves theta_ref_kp from parameters{} to data{} (named
// theta_ref_kp_data) and removes the prior on theta_ref_kp /
// mu_theta_ref_kp / sigma_theta_ref_kp from model{}, leaving the rest
// of the body bit-identical.

data {
  int<lower=1> n;
  int<lower=2> K;
  int<lower=2> p;

  // family_id_k_vector is homogeneous across slots in the initial
  // iteration of 8.6.D (heterogeneous slot families deferred). The
  // dispatcher uses family_id_k_vector[1] as the canonical family
  // identifier; the remaining slots inherit by family canonical
  // structure.
  array[K] int<lower=1, upper=13> family_id_k_vector;

  // inv_link_id_per_slot[k] is the canonical Stan inverse link of
  // slot k of the resolved family (inherited from Path B K dispatch).
  // The current iteration of 8.6.D only consumes inv_link_id_per_slot[2]
  // (log -> exp) for the dispersion / phi slot.
  array[K] int<lower=0, upper=2> inv_link_id_per_slot;

  // Per-slot use flags. The flat-pack offsets respect these flags;
  // when the slot is inactive, J_a_per_kp / J_b_per_kp entries are
  // zero and the dispatch loop reduces to a zero-cost branch.
  array[K] int<lower=0, upper=1> use_a_k;
  array[K] int<lower=0, upper=1> use_b_k;
  // use_W enforced to 0 upstream (D39 first iteration of 8.6.D).
  int<lower=0, upper=0> use_W;

  // Per-slot per-coord free-column counts.
  int<lower=0> J_a_max;
  int<lower=0> J_b_max;
  array[K, p] int<lower=0> J_a_per_kp;
  array[K, p] int<lower=0> J_b_per_kp;

  // Per-slot per-coord padded design matrices. Z_a_kp[k, j][i, 1:J_a_per_kp[k, j]]
  // and Z_b_kp[k, j][i, 1:J_b_per_kp[k, j]] are referenced; the
  // remaining padded columns are never touched.
  array[K, p] matrix[n, J_a_max] Z_a_kp;
  array[K, p] matrix[n, J_b_max] Z_b_kp;

  // Outcomes (matrix-column): y_real for continuous families
  // (Gaussian), y_int for count families (NB). Only the one matching
  // family_id_k_vector[1] is consumed.
  matrix[n, p] y_real;
  array[n, p] int y_int;

  // Per-slot per-coord anchor.
  array[K] vector[p] theta_anchor_kp;

  // Per-slot population-scope dispersion flags (derived on the R side
  // from family_id_k_vector and validated upstream).
  array[K] int<lower=0, upper=1> use_dispersion_y_k;
  array[K] int<lower=0, upper=1> use_dispersion_phi_k;

  // Grouping (Block 6.5, propagated to 8.6.D first iteration).
  int<lower=0, upper=1> use_groups;
  int<lower=1> J_groups;
  array[n] int<lower=1, upper=J_groups> group_id;

  // EB convention: every EB template declares K_slots and p_dim.
  // Path C: K_slots == K and p_dim == p (both strictly greater than 1).
  int<lower=2> K_slots;
  int<lower=2> p_dim;
}

transformed data {
  array[K, p] int J_a_free_kp;
  array[K, p] int J_b_free_kp;
  array[K, p] int a_raw_offset_kp;
  array[K, p] int c_b_raw_offset_kp;
  int total_J_a_free = 0;
  int total_J_b_free = 0;
  int any_use_a = 0;
  int any_use_b = 0;
  int sigma_y_size = 0;
  int phi_size = 0;
  array[K] int sigma_y_offset_k;
  array[K] int phi_offset_k;
  for (k in 1:K) {
    for (j in 1:p) {
      J_a_free_kp[k, j] = (use_a_k[k] == 1 && J_a_per_kp[k, j] > 0)
                            ? J_a_per_kp[k, j] : 0;
      J_b_free_kp[k, j] = (use_b_k[k] == 1 && J_b_per_kp[k, j] > 0)
                            ? J_b_per_kp[k, j] : 0;
      a_raw_offset_kp[k, j] = total_J_a_free;
      c_b_raw_offset_kp[k, j] = total_J_b_free;
      total_J_a_free += J_a_free_kp[k, j];
      total_J_b_free += J_b_free_kp[k, j];
      if (J_a_free_kp[k, j] > 0) any_use_a = 1;
      if (J_b_free_kp[k, j] > 0) any_use_b = 1;
    }
    sigma_y_offset_k[k] = sigma_y_size;
    phi_offset_k[k] = phi_size;
    sigma_y_size += use_dispersion_y_k[k];
    phi_size += use_dispersion_phi_k[k];
  }
  // Option A (RG.6, D96): compact sigma_a_k to slots carrying free a
  // coefficients (sum_j J_a_free_kp[k, j] > 0). A slot whose a()
  // component contributes no free coefficient --- no a() at all
  // (use_a_k == 0) or an intercept-only a() (a = ~ 1) absorbed into
  // theta_ref --- would otherwise declare a sampled-but-unused
  // half-prior scale: a flat direction (non-identification). sigma_a_idx
  // maps each slot to its position in the compacted vector (0 if none).
  // When every slot carries free a coefficients, n_sigma_a == K and
  // sigma_a_idx is the identity (the model stays bit-identical).
  int n_sigma_a = 0;
  array[K] int sigma_a_idx;
  for (k in 1:K) {
    int slot_free_a = 0;
    for (j in 1:p) {
      slot_free_a += J_a_free_kp[k, j];
    }
    if (slot_free_a > 0) {
      n_sigma_a += 1;
      sigma_a_idx[k] = n_sigma_a;
    } else {
      sigma_a_idx[k] = 0;
    }
  }
}

parameters {
  // Per-group, per-slot, per-coord anchor.
  array[J_groups, K] vector[p] theta_ref_kp;

  // Hyperparameters of the per-group anchor prior.
  array[use_groups, K] vector[p] mu_theta_ref_kp;
  array[use_groups, K] vector<lower=0>[p] sigma_theta_ref_kp;

  // Per-slot scale of the additive component (compacted to slots with
  // free a coefficients; see sigma_a_idx in transformed data, D96).
  array[n_sigma_a] real<lower=0> sigma_a_k;
  vector[total_J_a_free] a_raw;

  // Per-slot scale of the multiplicative component.
  array[K * any_use_b] real<lower=0> sigma_b_k;
  vector[total_J_b_free] c_b_kp_raw;

  // Per-slot dispersion entries.
  array[sigma_y_size] real<lower=0> sigma_y_pop_k;
  array[phi_size] real<lower=0> phi_pop_k;
}

transformed parameters {
  array[K] matrix[n, p] eta_kp;
  array[K, p] vector[J_a_max] a_coef_kp;
  array[K, p] vector[J_b_max] c_b_kp;
  array[K, p] vector[J_b_max] b_coef_kp;

  for (k in 1:K) {
    for (j in 1:p) {
      a_coef_kp[k, j] = rep_vector(0, J_a_max);
      c_b_kp[k, j] = rep_vector(0, J_b_max);
      b_coef_kp[k, j] = rep_vector(0, J_b_max);
    }
  }

  if (any_use_a == 1) {
    for (k in 1:K) {
      if (use_a_k[k] == 1) {
        for (j in 1:p) {
          if (J_a_free_kp[k, j] > 0) {
            for (jj in 1:J_a_free_kp[k, j]) {
              a_coef_kp[k, j][jj] = a_raw[a_raw_offset_kp[k, j] + jj]
                                      * sigma_a_k[sigma_a_idx[k]];
            }
          }
        }
      }
    }
  }

  if (any_use_b == 1) {
    for (k in 1:K) {
      if (use_b_k[k] == 1) {
        for (j in 1:p) {
          if (J_b_free_kp[k, j] > 0) {
            for (jj in 1:J_b_free_kp[k, j]) {
              c_b_kp[k, j][jj] = c_b_kp_raw[c_b_raw_offset_kp[k, j] + jj]
                                  * sigma_b_k[k];
            }
            // b_coef_kp[k, j] = c_b_kp[k, j] / theta_ref_kp[1, k][j]
            // is reported as a derived quantity only in the
            // single-anchor regime (use_groups == 0).
            if (use_groups == 0 && theta_ref_kp[1, k][j] != 0) {
              for (jj in 1:J_b_per_kp[k, j]) {
                b_coef_kp[k, j][jj] = c_b_kp[k, j][jj]
                                        / theta_ref_kp[1, k][j];
              }
            }
          }
        }
      }
    }
  }

  for (k in 1:K) {
    for (i in 1:n) {
      int g_i = group_id[i];
      for (j in 1:p) {
        real theta_ref_ikj = theta_ref_kp[g_i, k][j];
        real eta_ikj = theta_ref_ikj;
        if (use_a_k[k] == 1 && J_a_per_kp[k, j] > 0) {
          eta_ikj += dot_product(
            to_vector(Z_a_kp[k, j][i, 1:J_a_per_kp[k, j]]),
            a_coef_kp[k, j][1:J_a_per_kp[k, j]]
          );
        }
        if (use_b_k[k] == 1 && J_b_per_kp[k, j] > 0) {
          eta_ikj += dot_product(
            to_vector(Z_b_kp[k, j][i, 1:J_b_per_kp[k, j]]),
            c_b_kp[k, j][1:J_b_per_kp[k, j]]
          );
        }
        eta_kp[k][i, j] = eta_ikj;
      }
    }
  }
}

model {
  // Per-slot, per-coord anchor priors.
  if (use_groups == 1) {
    for (k in 1:K) {
      mu_theta_ref_kp[1, k] ~ {{PRIOR_THETA_REF}};
      sigma_theta_ref_kp[1, k] ~ {{PRIOR_SIGMA_THETA_REF}};
      for (g in 1:J_groups) {
        theta_ref_kp[g, k] ~ normal(mu_theta_ref_kp[1, k],
                                    sigma_theta_ref_kp[1, k]);
      }
    }
  } else {
    for (k in 1:K) {
      theta_ref_kp[1, k] ~ {{PRIOR_THETA_REF}};
    }
  }

  // Per-slot scale priors.
  if (any_use_a == 1) {
    sigma_a_k ~ {{PRIOR_SIGMA_A}};
    a_raw ~ normal(0, 1);
  }
  if (any_use_b == 1) {
    sigma_b_k ~ {{PRIOR_SIGMA_B}};
    c_b_kp_raw ~ normal(0, 1);
  }
  if (sigma_y_size > 0) {
    sigma_y_pop_k ~ {{PRIOR_SIGMA_Y}};
  }
  if (phi_size > 0) {
    phi_pop_k ~ {{PRIOR_PHI}};
  }

  // Likelihood: family dispatch + coord-wise factorization.
  // family_id_k_vector[1] selects the canonical Path B K=2 family;
  // slot 1 carries the location (mu under identity for Gaussian, log
  // for NB); slot 2 carries the dispersion (sigma for Gaussian via
  // inv_link 2 = exp, phi for NB via inv_link 2 = exp).
  if (family_id_k_vector[1] == 1) {
    // Gaussian K=2: y_real[i, j] ~ normal(eta_kp[1][i, j], sigma_ij)
    for (i in 1:n) {
      for (j in 1:p) {
        real sigma_ij = exp(eta_kp[2][i, j]);
        y_real[i, j] ~ normal(eta_kp[1][i, j], sigma_ij);
      }
    }
  } else if (family_id_k_vector[1] == 3) {
    // NB K=2: y_int[i, j] ~ neg_binomial_2(exp(eta_kp[1][i, j]), phi_ij)
    for (i in 1:n) {
      for (j in 1:p) {
        real phi_ij = exp(eta_kp[2][i, j]);
        y_int[i, j] ~ neg_binomial_2(exp(eta_kp[1][i, j]), phi_ij);
      }
    }
  }
  // family_id_k_vector[1] in {5, 6, 7, 8, 9, 10, 11, 12, 13}:
  // deferred to a later iteration of 8.6.D with the explicit
  // numerical caveat of opening Section 6.1. The dispatcher upstream
  // raises gdpar_unsupported_feature_error before reaching this
  // branch.
}

generated quantities {
  // Per-coord log-likelihood and posterior-predictive quantities.
  matrix[n, p] log_lik;
  matrix[n, p] theta_i;
  matrix[n, p] y_pred;

  if (family_id_k_vector[1] == 1) {
    for (i in 1:n) {
      for (j in 1:p) {
        real sigma_ij = exp(eta_kp[2][i, j]);
        theta_i[i, j] = eta_kp[1][i, j];
        log_lik[i, j] = normal_lpdf(y_real[i, j] | eta_kp[1][i, j],
                                     sigma_ij);
        y_pred[i, j] = normal_rng(eta_kp[1][i, j], sigma_ij);
      }
    }
  } else if (family_id_k_vector[1] == 3) {
    for (i in 1:n) {
      for (j in 1:p) {
        real phi_ij = exp(eta_kp[2][i, j]);
        real mu_ij = exp(eta_kp[1][i, j]);
        theta_i[i, j] = eta_kp[1][i, j];
        log_lik[i, j] = neg_binomial_2_lpmf(y_int[i, j] | mu_ij, phi_ij);
        y_pred[i, j] = neg_binomial_2_rng(mu_ij, phi_ij);
      }
    }
  }
}

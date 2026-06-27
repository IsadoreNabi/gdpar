// gdpar Path 1 -- Empirical Bayes (EB) Conditional Template, Path C
// (Sub-phase 8.6.D, decisions D36 = alpha + D37 = (i) canonized in
// Session 13a / 2026-05-25; decision D38'' = (h) canonized in
// Session 13b / 2026-05-25).
//
// See the marginal template amm_eb_marginal_KxP.stan header for the
// full architecture rationale (D38'' = (h) canonical Path C; Path A
// coord-wise factorization composed with Path B K-slot multi-
// parametric likelihood; coverage K=2 Gaussian / NB initial). This
// template moves theta_ref_kp from parameters{} to data{} (named
// theta_ref_kp_data) and removes its prior block from model{};
// otherwise the body is bit-identical to the marginal template.

data {
  int<lower=1> n;
  int<lower=2> K;
  int<lower=2> p;

  array[K] int<lower=1, upper=13> family_id_k_vector;
  array[K] int<lower=0, upper=2> inv_link_id_per_slot;

  array[K] int<lower=0, upper=1> use_a_k;
  array[K] int<lower=0, upper=1> use_b_k;
  int<lower=0, upper=0> use_W;

  int<lower=0> J_a_max;
  int<lower=0> J_b_max;
  array[K, p] int<lower=0> J_a_per_kp;
  array[K, p] int<lower=0> J_b_per_kp;

  array[K, p] matrix[n, J_a_max] Z_a_kp;
  array[K, p] matrix[n, J_b_max] Z_b_kp;

  matrix[n, p] y_real;
  array[n, p] int y_int;

  array[K] vector[p] theta_anchor_kp;

  array[K] int<lower=0, upper=1> use_dispersion_y_k;
  array[K] int<lower=0, upper=1> use_dispersion_phi_k;

  int<lower=0, upper=1> use_groups;
  int<lower=1> J_groups;
  array[n] int<lower=1, upper=J_groups> group_id;

  int<lower=2> K_slots;
  int<lower=2> p_dim;

  // Plug-in EB anchor passed via cmdstanr's automatic packing as an
  // R 3D array of shape [J_groups, K, p].
  array[J_groups, K] vector[p] theta_ref_kp_data;
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
  // theta_ref_kp is now plug-in EB data (theta_ref_kp_data).
  array[n_sigma_a] real<lower=0> sigma_a_k;
  vector[total_J_a_free] a_raw;

  array[K * any_use_b] real<lower=0> sigma_b_k;
  vector[total_J_b_free] c_b_kp_raw;

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
            if (use_groups == 0 && theta_ref_kp_data[1, k][j] != 0) {
              for (jj in 1:J_b_per_kp[k, j]) {
                b_coef_kp[k, j][jj] = c_b_kp[k, j][jj]
                                        / theta_ref_kp_data[1, k][j];
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
        real theta_ref_ikj = theta_ref_kp_data[g_i, k][j];
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
  // theta_ref_kp prior block intentionally absent: the anchor is data
  // here. The remaining priors mirror the marginal template.
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

  if (family_id_k_vector[1] == 1) {
    for (i in 1:n) {
      for (j in 1:p) {
        real sigma_ij = exp(eta_kp[2][i, j]);
        y_real[i, j] ~ normal(eta_kp[1][i, j], sigma_ij);
      }
    }
  } else if (family_id_k_vector[1] == 3) {
    for (i in 1:n) {
      for (j in 1:p) {
        real phi_ij = exp(eta_kp[2][i, j]);
        y_int[i, j] ~ neg_binomial_2(exp(eta_kp[1][i, j]), phi_ij);
      }
    }
  }
}

generated quantities {
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

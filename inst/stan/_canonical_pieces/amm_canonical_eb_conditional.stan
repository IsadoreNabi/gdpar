// gdpar Path 1 -- Empirical Bayes (EB) Conditional Template (Sub-phase 8.6.B)
//
// Used by R/eb.R via cmdstanr::cmdstan_model(...)$sample(...) in Step
// (iii) of the EB workflow described in v07 §11.1 (operational) and
// v07b §§4, 6 (theoretical):
//
//   Step (iii): plug the Step (i) EB point estimate
//               $\widehat\theta_{ref}^{EB}$ in as DATA and sample the
//               conditional posterior of xi = (a_coef, b_coef, W_raw,
//               sigma_*, phi) via HMC, returning $\Pi_n^{EB}(\xi \mid
//               \widehat\theta_{ref}^{EB})$ of v07 §2.1.
//
// Structural difference from amm_eb_marginal.stan (and amm_main.stan):
//   - theta_ref is moved from the parameters{} block to the data{}
//     block (here named theta_ref_data to make the regime explicit);
//   - the hyperparameters mu_theta_ref and sigma_theta_ref are no
//     longer sampled (they were exclusively a device to give theta_ref
//     a prior in the FB / EB-marginal regimes);
//   - the prior on theta_ref is removed from the model{} block (since
//     theta_ref is now fixed).
//
// Everything else is bit-identical to amm_eb_marginal.stan so that the
// conditional posterior of xi, given the plug-in
// $\widehat\theta_{ref}^{EB}$, matches the conditional posterior that
// the EB asymptotic theory of v07 §5 / v07b §4 analyses --- there is
// no extra modelling assumption between Steps (i) and (iii).
//
// Multivariate posture (Charter §2.7 of Sub-phase 8.6).
//   8.6.B operates under K = 1, p = 1, J >= 1 (degenerate scalar
//   single slot per group). The data block declares dimension fields
//   K_slots and p_dim with default 1; the body assumes K_slots == 1
//   and p_dim == 1 throughout. 8.6.C will extend the body by adding
//   K > 1 and p > 1 branches without rewriting the data declaration.

functions {
  // {{CANONICAL_HELPERS}}
}

data {
  int<lower=1> n;
  int<lower=1, upper=4> family_id;

  int<lower=0, upper=1> use_a;
  int<lower=0, upper=1> use_b;
  int<lower=0, upper=1> use_W;

  int<lower=0> J_a;
  int<lower=0> J_b;
  int<lower=0> dim_W;
  int<lower=0> d;

  matrix[n, J_a] Z_a;
  matrix[n, J_b] Z_b;
  matrix[n, d] X;

  vector[n] y_real;
  array[n] int y_int;

  real theta_anchor;

  int<lower=0, upper=1> use_dispersion_y;
  int<lower=0, upper=1> use_dispersion_phi;

  int<lower=0, upper=2> W_type_id;
  int<lower=0> W_n_knots_full;
  vector[W_n_knots_full] W_knots_full;
  int<lower=0> W_degree;

  int<lower=0, upper=1> use_groups;
  int<lower=1> J_groups;
  array[n] int<lower=1, upper=J_groups> group_id;

  // Sub-phase 8.6.B canonical-but-degenerate dimension fields. K_slots
  // and p_dim default to 1 for the 8.6.B regime and are checked R-side
  // by .gdpar_eb_validate_inputs(); 8.6.C will relax the K = 1 and
  // p = 1 guards and exercise these fields nontrivially.
  int<lower=1> K_slots;
  int<lower=1> p_dim;

  // EB conditional plug-in. The R-side helper
  // .gdpar_eb_maximize_marginal() populates theta_ref_data with
  // $\widehat\theta_{ref}^{EB}$ extracted from cmdstanr::laplace() on
  // amm_eb_marginal.stan; see R/eb.R for the projection and anti-
  // fragility logic of Charter §2.8.
  vector[J_groups] theta_ref_data;
}

transformed data {
  int J_a_free = (use_a == 1 && J_a > 0) ? J_a : 0;
  int J_b_free = (use_b == 1 && J_b > 0) ? J_b : 0;
  int dim_W_eff = (use_W == 1) ? dim_W : 0;
  int d_eff = (use_W == 1) ? d : 0;
}

parameters {
  // theta_ref is now data (theta_ref_data above); no hyperparameters
  // mu_theta_ref or sigma_theta_ref are sampled in the conditional
  // regime --- they were a device to give theta_ref a prior in the
  // marginal/FB regimes only.

  array[use_a == 1 && J_a > 0 ? 1 : 0] real<lower=0> sigma_a;
  vector[J_a_free] a_raw;

  array[use_b == 1 && J_b > 0 ? 1 : 0] real<lower=0> sigma_b;
  vector[J_b_free] c_b_raw;

  array[use_W == 1 && dim_W > 0 ? 1 : 0] real<lower=0> sigma_W;
  matrix[dim_W_eff, d_eff] W_raw;

  array[use_dispersion_y] real<lower=0> sigma_y;
  array[use_dispersion_phi] real<lower=0> phi;
}

transformed parameters {
  vector[J_a] a_coef = rep_vector(0, J_a);
  vector[J_b] c_b = rep_vector(0, J_b);
  vector[J_b] b_coef = rep_vector(0, J_b);
  vector[n] eta;

  if (use_a == 1 && J_a > 0) {
    for (j in 1:J_a_free) {
      a_coef[j] = a_raw[j]{{A_SCALE}};
    }
  }

  if (use_b == 1 && J_b > 0) {
    for (j in 1:J_b_free) {
      c_b[j] = c_b_raw[j] * sigma_b[1];
    }
    if (use_groups == 0 && theta_ref_data[1] != 0) {
      for (j in 1:J_b) {
        b_coef[j] = c_b[j] / theta_ref_data[1];
      }
    }
  }

  for (i in 1:n) {
    real theta_ref_i = theta_ref_data[group_id[i]];
    real eta_i = theta_ref_i;
    if (use_a == 1 && J_a > 0) {
      eta_i += Z_a[i] * a_coef;
    }
    if (use_b == 1 && J_b > 0) {
      eta_i += Z_b[i] * c_b;
    }
    if (use_W == 1 && dim_W > 0 && d > 0) {
      vector[d] W_diff_x = rep_vector(0, d);
      vector[dim_W] basis_diff = apply_W_basis_diff(
        W_type_id, theta_ref_i, theta_anchor, dim_W, W_degree,
        W_n_knots_full, W_knots_full
      );
      for (k in 1:dim_W) {
        for (jj in 1:d) {
          W_diff_x[jj] += basis_diff[k] * W_raw[k, jj]{{W_SCALE}};
        }
      }
      eta_i += dot_product(W_diff_x, to_vector(X[i]));
    }
    eta[i] = eta_i;
  }
}

model {
  // BEGIN PRIORS
  // No prior on theta_ref: theta_ref is data in the conditional
  // template. The hyperprior block of the marginal template
  // (mu_theta_ref, sigma_theta_ref) is removed here for the same
  // reason. xi components retain their canonical priors.
  if (use_a == 1 && J_a > 0) {
    sigma_a[1] ~ {{PRIOR_SIGMA_A}};
    a_raw ~ {{A_PRIOR}};
  }
  if (use_b == 1 && J_b > 0) {
    sigma_b[1] ~ {{PRIOR_SIGMA_B}};
    c_b_raw ~ normal(0, 1);
  }
  if (use_W == 1 && dim_W > 0) {
    sigma_W[1] ~ {{PRIOR_SIGMA_W}};
    to_vector(W_raw) ~ {{W_PRIOR}};
  }
  if (use_dispersion_y == 1) {
    sigma_y[1] ~ {{PRIOR_SIGMA_Y}};
  }
  if (use_dispersion_phi == 1) {
    phi[1] ~ {{PRIOR_PHI}};
  }
  // END PRIORS

  if (family_id == 1) {
    y_real ~ normal(eta, sigma_y[1]);
  } else if (family_id == 2) {
    y_int ~ poisson_log(eta);
  } else if (family_id == 3) {
    y_int ~ neg_binomial_2_log(eta, phi[1]);
  } else if (family_id == 4) {
    y_int ~ bernoulli_logit(eta);
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] theta_i = eta;
  vector[n] y_pred;

  for (i in 1:n) {
    if (family_id == 1) {
      log_lik[i] = normal_lpdf(y_real[i] | eta[i], sigma_y[1]);
      y_pred[i] = normal_rng(eta[i], sigma_y[1]);
    } else if (family_id == 2) {
      log_lik[i] = poisson_log_lpmf(y_int[i] | eta[i]);
      y_pred[i] = poisson_log_rng(eta[i]);
    } else if (family_id == 3) {
      log_lik[i] = neg_binomial_2_log_lpmf(y_int[i] | eta[i], phi[1]);
      y_pred[i] = neg_binomial_2_log_rng(eta[i], phi[1]);
    } else if (family_id == 4) {
      log_lik[i] = bernoulli_logit_lpmf(y_int[i] | eta[i]);
      y_pred[i] = bernoulli_logit_rng(eta[i]);
    }
  }
}

// Sub-phase 8.3.8 (2026-05-22): Cox-de Boor B-spline basis evaluation
  // for the W modulating component. The augmented knot vector
  // W_knots_full (length n_int + 2*(W_degree+1)) is constructed R-side
  // by .gdpar_bspline_knots_full() and shipped as data; the recursion
  // re-evaluates the basis at theta_ref_i and theta_anchor in every
  // Hamiltonian Monte Carlo step (these are parameters, so a R-side
  // pre-computed matrix would not be differentiable in the sampler).
  // The returned dimension W_per_k_dim equals the total number of
  // Cox-de Boor B-splines (n_int + W_degree + 1) minus one, mirroring
  // the dimensional convention of splines::bs(intercept = FALSE); the
  // first basis is dropped. The exact ordering does not need to bit-
  // match splines::bs because W_raw absorbs any orthogonal rotation
  // under sampling, but the dimension count must agree (R-side dim is
  // computed by materialize_W_basis() as length(knots) + degree per
  // coordinate; for the K = 1 + p = 1 template W_per_k_dim is dim_W).
  vector bspline_basis_eval(real x, int W_per_k_dim, int W_degree,
                            int W_n_knots_full, vector W_knots_full) {
    int n_basis_total = W_n_knots_full - W_degree - 1;
    int n_iter_max = n_basis_total + W_degree;
    vector[n_iter_max] B;
    for (i in 1:n_iter_max) {
      if (W_knots_full[i] == W_knots_full[i + 1]) {
        B[i] = 0.0;
      } else if (i == n_iter_max && x >= W_knots_full[i + 1]) {
        B[i] = 1.0;
      } else {
        B[i] = (x >= W_knots_full[i] && x < W_knots_full[i + 1]) ? 1.0 : 0.0;
      }
    }
    for (d_ord in 1:W_degree) {
      int n_iter = n_basis_total + W_degree - d_ord;
      vector[n_iter] B_new;
      for (i in 1:n_iter) {
        real denom1 = W_knots_full[i + d_ord] - W_knots_full[i];
        real denom2 = W_knots_full[i + d_ord + 1] - W_knots_full[i + 1];
        real term1 = (denom1 > 0.0)
          ? (x - W_knots_full[i]) / denom1 * B[i] : 0.0;
        real term2 = (denom2 > 0.0)
          ? (W_knots_full[i + d_ord + 1] - x) / denom2 * B[i + 1] : 0.0;
        B_new[i] = term1 + term2;
      }
      for (i in 1:n_iter) B[i] = B_new[i];
    }
    vector[W_per_k_dim] out;
    for (j in 1:W_per_k_dim) out[j] = B[j + 1];
    return out;
  }

  // Sub-phase 8.3.8: dispatcher producing the anchored basis difference
  // basis(theta) - basis(anchor) per W_type_id. Returns rep_vector(0,
  // W_per_k_dim) when basis_type_id == 0 (W off, defensive).
  vector apply_W_basis_diff(int basis_type_id, real theta, real anchor,
                            int W_per_k_dim, int W_degree,
                            int W_n_knots_full, vector W_knots_full) {
    vector[W_per_k_dim] diff = rep_vector(0.0, W_per_k_dim);
    if (basis_type_id == 1) {
      for (jj in 1:W_per_k_dim) {
        diff[jj] = pow(theta, jj) - pow(anchor, jj);
      }
    } else if (basis_type_id == 2) {
      diff = bspline_basis_eval(theta, W_per_k_dim, W_degree,
                                W_n_knots_full, W_knots_full)
           - bspline_basis_eval(anchor, W_per_k_dim, W_degree,
                                W_n_knots_full, W_knots_full);
    }
    return diff;
  }
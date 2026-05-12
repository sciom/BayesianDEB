// ============================================================
// BayesianDEB: DEBtox model — TKTD with DEB energetics
// 4-state: reserve E, structure V, repro buffer R, scaled damage Dw
// Stress on assimilation (mode of action)
// Within-chain parallelism via reduce_sum
// Supports Arrhenius temperature correction (optional)
// Rate parameters corrected: p_Am, p_M, v, k_d
// NOT corrected: kappa, E_G, z_w, b_w
// ============================================================

functions {
  vector deb_tox_ode(real t, vector x,
                     real p_Am, real p_M, real kappa,
                     real v_, real E_G,
                     real k_d, real z_w, real b_w,
                     real f, real C_w) {
    real E  = x[1];
    real V  = x[2];
    real R  = x[3];
    real Dw = x[4];

    real dDw = k_d * (fmax(C_w - z_w, 0.0) - Dw);
    real s = b_w * fmax(Dw, 0.0);

    real L = pow(V, 1.0 / 3.0);

    real p_A = f * p_Am * L * L * fmax(1.0 - s, 0.0);
    real p_C = E * v_ * L / (E + E_G * V + 1e-12);
    real p_M_flux = p_M * V;

    real dE = p_A - p_C;
    real dV = (kappa * p_C - p_M_flux) / E_G;
    real dR = fmax((1 - kappa) * p_C, 0.0);

    if (V < 1e-12 && dV < 0) dV = 0;

    return [dE, dV, dR, dDw]';
  }

  real partial_log_lik(array[] int grp_slice,
                       int start, int end,
                       array[] real C_w,
                       int max_N_obs,
                       array[] int N_obs,
                       array[,] real t_obs,
                       array[,] real L_obs,
                       real f_food,
                       int has_repro,
                       array[] int N_R,
                       int max_N_R,
                       array[,] int R_counts,
                       array[,] int idx_R_start,
                       array[,] int idx_R_end,
                       real p_Am, real p_M, real kappa,
                       real v, real E_G,
                       real E0, real L0, real sigma_L,
                       real k_d, real z_w, real b_w,
                       real k_R, real phi_R,
                       real c_T,
                       int obs_growth, real obs_nu,
                       int obs_repro) {
    real ll = 0;
    for (k in 1:size(grp_slice)) {
      int g = grp_slice[k];

      real V0 = L0 * L0 * L0;
      vector[4] x0 = [E0 * V0, V0, 0.0, 0.0]';

      array[N_obs[g]] real t_g;
      for (i in 1:N_obs[g]) t_g[i] = t_obs[g, i];

      array[N_obs[g]] vector[4] x_sol = ode_bdf_tol(
        deb_tox_ode, x0, 0.0, t_g,
        1e-5, 1e-8, 1000000,
        p_Am * c_T, p_M * c_T, kappa, v * c_T, E_G,
        k_d * c_T, z_w, b_w,
        f_food, C_w[g]
      );

      // Growth likelihood
      for (i in 1:N_obs[g]) {
        real L_hat = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
        if (!is_nan(L_obs[g, i])) {
          if (obs_growth == 1)
            ll += normal_lpdf(L_obs[g, i] | L_hat, sigma_L);
          else if (obs_growth == 2)
            ll += lognormal_lpdf(L_obs[g, i] | log(fmax(L_hat, 1e-12)), sigma_L);
          else
            ll += student_t_lpdf(L_obs[g, i] | obs_nu, L_hat, sigma_L);
        }
      }

      // Reproduction likelihood (interval counts, consistent with growth_repro)
      if (has_repro == 1) {
        for (i in 1:N_R[g]) {
          real delta_R = k_R * fmax(
            x_sol[idx_R_end[g, i]][3] - x_sol[idx_R_start[g, i]][3], 0.0);
          if (obs_repro == 1)
            ll += neg_binomial_2_lpmf(R_counts[g, i] | fmax(delta_R, 1e-6), phi_R);
          else
            ll += poisson_lpmf(R_counts[g, i] | fmax(delta_R, 1e-6));
        }
      }
    }
    return ll;
  }
}

data {
  int<lower=1> N_groups;
  array[N_groups] real<lower=0> C_w;
  int<lower=1> max_N_obs;
  array[N_groups] int<lower=1> N_obs;
  array[N_groups, max_N_obs] real t_obs;
  array[N_groups, max_N_obs] real L_obs;
  real<lower=0, upper=1> f_food;

  int<lower=0> has_repro;
  array[N_groups] int<lower=0> N_R;
  int<lower=0> max_N_R;
  array[N_groups, max_N_R] int<lower=0> R_counts;
  array[N_groups, max_N_R] int<lower=1> idx_R_start;
  array[N_groups, max_N_R] int<lower=1> idx_R_end;

  // Observation models
  int<lower=1, upper=3> obs_growth;  // 1=normal, 2=lognormal, 3=student_t
  real<lower=1> obs_nu;
  int<lower=1, upper=2> obs_repro;   // 1=negbinom, 2=poisson

  // Temperature correction (optional)
  int<lower=0, upper=1> has_temperature;
  real<lower=0> T_obs;
  real<lower=0> T_ref;
  real<lower=0> T_A;

  // Prior hyperparameters
  real prior_p_Am_mu;
  real<lower=0> prior_p_Am_sd;
  real prior_p_M_mu;
  real<lower=0> prior_p_M_sd;
  real<lower=0> prior_kappa_a;
  real<lower=0> prior_kappa_b;
  real prior_v_mu;
  real<lower=0> prior_v_sd;
  real prior_E_G_mu;
  real<lower=0> prior_E_G_sd;
  real prior_E0_mu;
  real<lower=0> prior_E0_sd;
  real prior_L0_mu;
  real<lower=0> prior_L0_sd;
  real<lower=0> prior_sigma_L_sd;
  real prior_k_d_mu;
  real<lower=0> prior_k_d_sd;
  real prior_z_w_mu;
  real<lower=0> prior_z_w_sd;
  real prior_b_w_mu;
  real<lower=0> prior_b_w_sd;
  real prior_k_R_mu;
  real<lower=0> prior_k_R_sd;
  real prior_phi_R_mu;
  real<lower=0> prior_phi_R_sd;
}

transformed data {
  array[N_groups] int grp_idx;
  for (g in 1:N_groups) grp_idx[g] = g;

  real c_T = has_temperature
    ? exp(T_A / T_ref - T_A / T_obs)
    : 1.0;
}

parameters {
  real<lower=0> p_Am;
  real<lower=0> p_M;
  real<lower=0, upper=1> kappa;
  real<lower=0> v;
  real<lower=0> E_G;
  real<lower=0> E0;
  real<lower=0> L0;
  real<lower=0> sigma_L;

  real<lower=0> k_d;
  real<lower=0> z_w;
  real<lower=0> b_w;
  real<lower=0> k_R;    // repro buffer to offspring conversion
  real<lower=0> phi_R;  // NegBin overdispersion (only used if obs_repro==1)
}

model {
  p_Am  ~ lognormal(prior_p_Am_mu, prior_p_Am_sd);
  p_M   ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa ~ beta(prior_kappa_a, prior_kappa_b);
  v     ~ lognormal(prior_v_mu, prior_v_sd);
  E_G   ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  E0    ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0    ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);  // half-normal via <lower=0> constraint

  k_d ~ lognormal(prior_k_d_mu, prior_k_d_sd);
  z_w ~ lognormal(prior_z_w_mu, prior_z_w_sd);
  b_w ~ lognormal(prior_b_w_mu, prior_b_w_sd);
  k_R   ~ lognormal(prior_k_R_mu, prior_k_R_sd);
  phi_R ~ lognormal(prior_phi_R_mu, prior_phi_R_sd);

  target += reduce_sum(partial_log_lik, grp_idx, 1,
                       C_w, max_N_obs, N_obs, t_obs, L_obs, f_food,
                       has_repro, N_R, max_N_R, R_counts,
                       idx_R_start, idx_R_end,
                       p_Am, p_M, kappa, v, E_G,
                       E0, L0, sigma_L,
                       k_d, z_w, b_w, k_R, phi_R,
                       c_T, obs_growth, obs_nu, obs_repro);
}

generated quantities {
  real EC50 = z_w + 0.5 / b_w;
  real NEC  = z_w;
}

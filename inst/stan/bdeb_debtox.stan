// ============================================================
// BayesianDEB: DEBtox model — TKTD with DEB energetics
// 4-state: reserve E, structure V, repro buffer R, scaled damage Dw
// Stress on assimilation (mode of action)
// ============================================================

functions {
  vector deb_tox_ode(real t, vector x, array[] real theta,
                     array[] real env) {
    real E  = x[1];
    real V  = x[2];
    real R  = x[3];
    real Dw = x[4];

    real p_Am  = theta[1];
    real p_M   = theta[2];
    real kappa = theta[3];
    real v_    = theta[4];
    real E_G   = theta[5];
    real k_d   = theta[6];
    real z_w   = theta[7];
    real b_w   = theta[8];

    real f   = env[1];
    real C_w = env[2];

    // Toxicokinetics: scaled damage dynamics
    real dDw = k_d * (fmax(C_w - z_w, 0.0) - Dw);
    real s = b_w * fmax(Dw, 0.0);

    real L = pow(V, 1.0 / 3.0);

    // Stress on assimilation
    real p_A = f * p_Am * L * L * fmax(1.0 - s, 0.0);
    real p_C = E * v_ * L / (E + E_G * V + 1e-12);
    real p_M_flux = p_M * V;

    real dE = p_A - p_C;
    real dV = (kappa * p_C - p_M_flux) / E_G;
    real dR = fmax((1 - kappa) * p_C, 0.0);

    if (V < 1e-12 && dV < 0) dV = 0;

    return [dE, dV, dR, dDw]';
  }
}

data {
  int<lower=1> N_groups;                   // number of concentration groups
  array[N_groups] real<lower=0> C_w;       // external concentrations
  int<lower=1> max_N_obs;
  array[N_groups] int<lower=1> N_obs;
  array[N_groups, max_N_obs] real t_obs;
  array[N_groups, max_N_obs] real L_obs;   // NaN for missing
  real<lower=0, upper=1> f_food;

  // Reproduction data (optional, per group)
  int<lower=0> has_repro;
  array[N_groups] int<lower=0> N_R;
  int<lower=0> max_N_R;
  array[N_groups, max_N_R] int<lower=0> R_counts;
  array[N_groups, max_N_R] int<lower=1> idx_R_end;

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

  // Tox parameters
  real<lower=0> k_d;
  real<lower=0> z_w;       // NEC (no-effect concentration)
  real<lower=0> b_w;       // effect intensity
}

model {
  // Priors — DEB
  p_Am  ~ lognormal(prior_p_Am_mu, prior_p_Am_sd);
  p_M   ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa ~ beta(prior_kappa_a, prior_kappa_b);
  v     ~ lognormal(prior_v_mu, prior_v_sd);
  E_G   ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  E0    ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0    ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);

  // Priors — Tox
  k_d ~ lognormal(prior_k_d_mu, prior_k_d_sd);
  z_w ~ lognormal(prior_z_w_mu, prior_z_w_sd);
  b_w ~ lognormal(prior_b_w_mu, prior_b_w_sd);

  // Per-group ODE solve and likelihood
  for (g in 1:N_groups) {
    vector[4] x0;
    real V0 = L0 * L0 * L0;
    x0[1] = E0 * V0;
    x0[2] = V0;
    x0[3] = 0.0;
    x0[4] = 0.0;

    array[8] real theta = {p_Am, p_M, kappa, v, E_G, k_d, z_w, b_w};
    array[2] real env = {f_food, C_w[g]};

    array[N_obs[g]] real t_g;
    for (i in 1:N_obs[g]) t_g[i] = t_obs[g, i];

    array[N_obs[g]] vector[4] x_sol = ode_bdf(
      deb_tox_ode, x0, 0.0, t_g, theta, env, 1e-6, 1e-6, 1e4
    );

    // Growth likelihood
    for (i in 1:N_obs[g]) {
      real L_hat = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
      if (!is_nan(L_obs[g, i])) {
        L_obs[g, i] ~ normal(L_hat, sigma_L);
      }
    }

    // Reproduction likelihood (if present)
    if (has_repro == 1) {
      for (i in 1:N_R[g]) {
        real R_pred = fmax(x_sol[idx_R_end[g, i]][3], 1e-6);
        R_counts[g, i] ~ poisson(R_pred);
      }
    }
  }
}

generated quantities {
  // EC50 approximation: concentration where stress s = 0.5
  // s = b_w * max(C - z_w, 0) at steady state (Dw = C - z_w)
  // s = 0.5 => C = z_w + 0.5 / b_w
  real EC50 = z_w + 0.5 / b_w;
  real NEC  = z_w;
}

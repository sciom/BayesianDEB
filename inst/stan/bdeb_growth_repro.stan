// ============================================================
// BayesianDEB: Individual growth + reproduction model
// 3-state DEB (E, V, R) with Gaussian length + NegBin reproduction
// ============================================================

functions {
  vector deb_growth_repro_ode(real t, vector x, array[] real theta,
                              array[] real env) {
    real E  = x[1];
    real V  = x[2];
    real R  = x[3];

    real p_Am  = theta[1];
    real p_M   = theta[2];
    real kappa = theta[3];
    real v_    = theta[4];
    real E_G   = theta[5];
    real k_J   = theta[6];
    real E_Hp  = theta[7];
    real f     = env[1];

    real L = pow(V, 1.0 / 3.0);
    real p_A = f * p_Am * L * L;
    real p_C = E * v_ * L / (E + E_G * V + 1e-12);
    real p_M_flux = p_M * V;
    real p_J = k_J * E_Hp;

    real dE = p_A - p_C;
    real dV = (kappa * p_C - p_M_flux) / E_G;
    real dR = fmax((1 - kappa) * p_C - p_J, 0.0);

    if (V < 1e-12 && dV < 0) dV = 0;

    return [dE, dV, dR]';
  }
}

data {
  // Growth data
  int<lower=0> N_L;
  array[N_L] real<lower=0> t_L;
  array[N_L] real<lower=0> L_obs;

  // Reproduction data (interval counts)
  int<lower=0> N_R;
  array[N_R] real<lower=0> t_R_start;
  array[N_R] real<lower=0> t_R_end;
  array[N_R] int<lower=0> R_counts;

  real<lower=0, upper=1> f_food;

  // Combined unique times for ODE solving
  int<lower=1> N_times;
  array[N_times] real<lower=0> t_all;
  // Index mapping: which t_all entry corresponds to each observation
  array[N_L] int<lower=1> idx_L;
  array[N_R] int<lower=1> idx_R_start;
  array[N_R] int<lower=1> idx_R_end;

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
  real prior_k_J_mu;
  real<lower=0> prior_k_J_sd;
  real prior_E_Hp_mu;
  real<lower=0> prior_E_Hp_sd;
  real prior_E0_mu;
  real<lower=0> prior_E0_sd;
  real prior_L0_mu;
  real<lower=0> prior_L0_sd;
  real<lower=0> prior_sigma_L_sd;
  real<lower=0> prior_k_R_mu;
  real<lower=0> prior_k_R_sd;
  real<lower=0> prior_phi_R_mu;
  real<lower=0> prior_phi_R_sd;
}

parameters {
  real<lower=0> p_Am;
  real<lower=0> p_M;
  real<lower=0, upper=1> kappa;
  real<lower=0> v;
  real<lower=0> E_G;
  real<lower=0> k_J;
  real<lower=0> E_Hp;
  real<lower=0> E0;
  real<lower=0> L0;
  real<lower=0> sigma_L;
  real<lower=0> k_R;       // repro buffer to offspring conversion
  real<lower=0> phi_R;     // NegBin overdispersion
}

transformed parameters {
  array[N_times] vector[3] x_sol;
  array[N_times] real L_hat;
  array[N_times] real R_hat;

  {
    vector[3] x0;
    real V0 = L0 * L0 * L0;
    x0[1] = E0 * V0;
    x0[2] = V0;
    x0[3] = 0.0;

    array[7] real theta = {p_Am, p_M, kappa, v, E_G, k_J, E_Hp};
    array[1] real env = {f_food};

    x_sol = ode_bdf(deb_growth_repro_ode, x0, 0.0, t_all, theta, env,
                    1e-6, 1e-6, 1e4);
  }

  for (i in 1:N_times) {
    L_hat[i] = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
    R_hat[i] = fmax(x_sol[i][3], 0.0);
  }
}

model {
  // Priors
  p_Am  ~ lognormal(prior_p_Am_mu, prior_p_Am_sd);
  p_M   ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa ~ beta(prior_kappa_a, prior_kappa_b);
  v     ~ lognormal(prior_v_mu, prior_v_sd);
  E_G   ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  k_J   ~ lognormal(prior_k_J_mu, prior_k_J_sd);
  E_Hp  ~ lognormal(prior_E_Hp_mu, prior_E_Hp_sd);
  E0    ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0    ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);
  k_R   ~ lognormal(prior_k_R_mu, prior_k_R_sd);
  phi_R ~ lognormal(prior_phi_R_mu, prior_phi_R_sd);

  // Growth likelihood
  for (i in 1:N_L) {
    L_obs[i] ~ normal(L_hat[idx_L[i]], sigma_L);
  }

  // Reproduction likelihood (interval counts, negative binomial)
  for (i in 1:N_R) {
    real delta_R = k_R * fmax(R_hat[idx_R_end[i]] - R_hat[idx_R_start[i]], 0.0);
    R_counts[i] ~ neg_binomial_2(fmax(delta_R, 1e-6), phi_R);
  }
}

generated quantities {
  array[N_L] real L_rep;
  array[N_R] int R_rep;
  array[N_L] real log_lik_L;
  array[N_R] real log_lik_R;

  for (i in 1:N_L) {
    L_rep[i] = normal_rng(L_hat[idx_L[i]], sigma_L);
    log_lik_L[i] = normal_lpdf(L_obs[i] | L_hat[idx_L[i]], sigma_L);
  }

  for (i in 1:N_R) {
    real delta_R = k_R * fmax(R_hat[idx_R_end[i]] - R_hat[idx_R_start[i]], 0.0);
    R_rep[i] = neg_binomial_2_rng(fmax(delta_R, 1e-6), phi_R);
    log_lik_R[i] = neg_binomial_2_lpmf(R_counts[i] | fmax(delta_R, 1e-6), phi_R);
  }
}

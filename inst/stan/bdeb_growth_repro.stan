// ============================================================
// BayesianDEB: Individual growth + reproduction model
// 3-state DEB (E, V, R)
// Growth obs: 1=normal, 2=lognormal, 3=student_t
// Repro obs:  1=negbinom, 2=poisson
// Supports Arrhenius temperature correction (optional)
// ============================================================

functions {
  vector deb_growth_repro_ode(real t, vector x,
                              real p_Am, real p_M, real kappa,
                              real v_, real E_G, real k_J,
                              real E_Hp, real f) {
    real E  = x[1];
    real V  = x[2];
    real R  = x[3];

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
  int<lower=0> N_L;
  array[N_L] real<lower=0> t_L;
  array[N_L] real<lower=0> L_obs;

  int<lower=0> N_R;
  array[N_R] real<lower=0> t_R_start;
  array[N_R] real<lower=0> t_R_end;
  array[N_R] int<lower=0> R_counts;

  real<lower=0, upper=1> f_food;

  int<lower=1> N_times;
  array[N_times] real<lower=0> t_all;
  array[N_L] int<lower=1> idx_L;
  array[N_R] int<lower=1> idx_R_start;
  array[N_R] int<lower=1> idx_R_end;

  // Observation models
  int<lower=1, upper=3> obs_growth;  // 1=normal, 2=lognormal, 3=student_t
  real<lower=1> obs_nu;
  int<lower=1, upper=2> obs_repro;   // 1=negbinom, 2=poisson

  // Temperature
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
  real prior_k_J_mu;
  real<lower=0> prior_k_J_sd;
  real prior_E_Hp_mu;
  real<lower=0> prior_E_Hp_sd;
  real prior_E0_mu;
  real<lower=0> prior_E0_sd;
  real prior_L0_mu;
  real<lower=0> prior_L0_sd;
  real<lower=0> prior_sigma_L_sd;
  real prior_k_R_mu;
  real<lower=0> prior_k_R_sd;
  real prior_phi_R_mu;
  real<lower=0> prior_phi_R_sd;
}

transformed data {
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
  real<lower=0> k_J;
  real<lower=0> E_Hp;
  real<lower=0> E0;
  real<lower=0> L0;
  real<lower=0> sigma_L;
  real<lower=0> k_R;
  real<lower=0> phi_R;     // only used if obs_repro == 1
}

transformed parameters {
  array[N_times] vector[3] x_sol;
  array[N_times] real L_hat;
  array[N_times] real R_hat;

  {
    real V0 = L0 * L0 * L0;
    vector[3] x0 = [E0 * V0, V0, 0.0]';

    x_sol = ode_bdf_tol(deb_growth_repro_ode, x0, 0.0, t_all,
                        1e-5, 1e-8, 1000000,
                        p_Am * c_T, p_M * c_T, kappa,
                        v * c_T, E_G, k_J * c_T, E_Hp, f_food);
  }

  for (i in 1:N_times) {
    L_hat[i] = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
    R_hat[i] = fmax(x_sol[i][3], 0.0);
  }
}

model {
  p_Am  ~ lognormal(prior_p_Am_mu, prior_p_Am_sd);
  p_M   ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa ~ beta(prior_kappa_a, prior_kappa_b);
  v     ~ lognormal(prior_v_mu, prior_v_sd);
  E_G   ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  k_J   ~ lognormal(prior_k_J_mu, prior_k_J_sd);
  E_Hp  ~ lognormal(prior_E_Hp_mu, prior_E_Hp_sd);
  E0    ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0    ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);  // half-normal via <lower=0> constraint
  k_R   ~ lognormal(prior_k_R_mu, prior_k_R_sd);
  phi_R ~ lognormal(prior_phi_R_mu, prior_phi_R_sd);

  // Growth likelihood
  for (i in 1:N_L) {
    if (obs_growth == 1)
      L_obs[i] ~ normal(L_hat[idx_L[i]], sigma_L);
    else if (obs_growth == 2)
      L_obs[i] ~ lognormal(log(fmax(L_hat[idx_L[i]], 1e-12)), sigma_L);
    else
      L_obs[i] ~ student_t(obs_nu, L_hat[idx_L[i]], sigma_L);
  }

  // Reproduction likelihood
  for (i in 1:N_R) {
    real delta_R = k_R * fmax(R_hat[idx_R_end[i]] - R_hat[idx_R_start[i]], 0.0);
    if (obs_repro == 1)
      R_counts[i] ~ neg_binomial_2(fmax(delta_R, 1e-6), phi_R);
    else
      R_counts[i] ~ poisson(fmax(delta_R, 1e-6));
  }
}

generated quantities {
  array[N_L] real L_rep;
  array[N_R] int R_rep;
  array[N_L] real log_lik_L;
  array[N_R] real log_lik_R;

  for (i in 1:N_L) {
    real mu_L = L_hat[idx_L[i]];
    if (obs_growth == 1) {
      L_rep[i] = normal_rng(mu_L, sigma_L);
      log_lik_L[i] = normal_lpdf(L_obs[i] | mu_L, sigma_L);
    } else if (obs_growth == 2) {
      L_rep[i] = lognormal_rng(log(fmax(mu_L, 1e-12)), sigma_L);
      log_lik_L[i] = lognormal_lpdf(L_obs[i] | log(fmax(mu_L, 1e-12)), sigma_L);
    } else {
      L_rep[i] = student_t_rng(obs_nu, mu_L, sigma_L);
      log_lik_L[i] = student_t_lpdf(L_obs[i] | obs_nu, mu_L, sigma_L);
    }
  }

  for (i in 1:N_R) {
    real delta_R = k_R * fmax(R_hat[idx_R_end[i]] - R_hat[idx_R_start[i]], 0.0);
    if (obs_repro == 1) {
      R_rep[i] = neg_binomial_2_rng(fmax(delta_R, 1e-6), phi_R);
      log_lik_R[i] = neg_binomial_2_lpmf(R_counts[i] | fmax(delta_R, 1e-6), phi_R);
    } else {
      R_rep[i] = poisson_rng(fmax(delta_R, 1e-6));
      log_lik_R[i] = poisson_lpmf(R_counts[i] | fmax(delta_R, 1e-6));
    }
  }
}

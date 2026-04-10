// ============================================================
// BayesianDEB: Individual-level growth model
// 2-state DEB (reserve E, structure V)
// Observation models: 1=normal, 2=lognormal, 3=student_t
// Supports Arrhenius temperature correction (optional)
// ============================================================

functions {
  vector deb_growth_ode(real t, vector x,
                        real p_Am, real p_M, real kappa,
                        real v_, real E_G, real f) {
    real E = x[1];
    real V = x[2];

    real L = pow(V, 1.0 / 3.0);
    real p_A = f * p_Am * L * L;
    real p_C = E * v_ * L / (E + E_G * V + 1e-12);
    real p_M_flux = p_M * V;

    real dE = p_A - p_C;
    real dV = (kappa * p_C - p_M_flux) / E_G;
    if (V < 1e-12 && dV < 0) dV = 0;

    return [dE, dV]';
  }
}

data {
  int<lower=1> N_obs;
  array[N_obs] real<lower=0> t_obs;
  array[N_obs] real<lower=0> L_obs;
  real<lower=0, upper=1> f_food;

  // Observation model: 1=normal, 2=lognormal, 3=student_t
  int<lower=1, upper=3> obs_growth;  // 1=normal, 2=lognormal, 3=student_t
  real<lower=1> obs_nu;              // df for student_t (only if obs_growth=3)
  int<lower=1, upper=2> obs_repro;   // unused here, declared for data compatibility

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
  real<lower=0> E0;
  real<lower=0> L0;
  real<lower=0> sigma_L;
}

transformed parameters {
  array[N_obs] vector[2] x_sol;
  array[N_obs] real L_hat;

  {
    real V0 = L0 * L0 * L0;
    vector[2] x0 = [E0 * V0, V0]';

    x_sol = ode_bdf_tol(deb_growth_ode, x0, 0.0, t_obs,
                        1e-6, 1e-6, 10000,
                        p_Am * c_T, p_M * c_T, kappa,
                        v * c_T, E_G, f_food);
  }

  for (i in 1:N_obs) {
    L_hat[i] = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
  }
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

  // Growth likelihood — switched by obs_growth flag
  if (obs_growth == 1) {
    L_obs ~ normal(L_hat, sigma_L);
  } else if (obs_growth == 2) {
    for (i in 1:N_obs)
      L_obs[i] ~ lognormal(log(fmax(L_hat[i], 1e-12)), sigma_L);
  } else {
    for (i in 1:N_obs)
      L_obs[i] ~ student_t(obs_nu, L_hat[i], sigma_L);
  }
}

generated quantities {
  array[N_obs] real L_rep;
  array[N_obs] real log_lik;

  for (i in 1:N_obs) {
    if (obs_growth == 1) {
      L_rep[i] = normal_rng(L_hat[i], sigma_L);
      log_lik[i] = normal_lpdf(L_obs[i] | L_hat[i], sigma_L);
    } else if (obs_growth == 2) {
      L_rep[i] = lognormal_rng(log(fmax(L_hat[i], 1e-12)), sigma_L);
      log_lik[i] = lognormal_lpdf(L_obs[i] | log(fmax(L_hat[i], 1e-12)), sigma_L);
    } else {
      L_rep[i] = student_t_rng(obs_nu, L_hat[i], sigma_L);
      log_lik[i] = student_t_lpdf(L_obs[i] | obs_nu, L_hat[i], sigma_L);
    }
  }
}

// ============================================================
// BayesianDEB: Individual-level growth model
// 2-state DEB (reserve E, structure V) with Gaussian observations
// ============================================================

functions {
  vector deb_growth_ode(real t, vector x, array[] real theta,
                        array[] real env) {
    real E = x[1];
    real V = x[2];
    real p_Am  = theta[1];
    real p_M   = theta[2];
    real kappa = theta[3];
    real v_    = theta[4];
    real E_G   = theta[5];
    real f     = env[1];

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
  real<lower=0, upper=1> f_food;        // scaled functional response

  // Prior hyperparameters (passed from R)
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
    vector[2] x0;
    real V0 = L0 * L0 * L0;
    x0[1] = E0 * V0;
    x0[2] = V0;

    array[5] real theta = {p_Am, p_M, kappa, v, E_G};
    array[1] real env = {f_food};

    x_sol = ode_bdf(deb_growth_ode, x0, 0.0, t_obs, theta, env,
                    1e-6, 1e-6, 1e4);
  }

  for (i in 1:N_obs) {
    L_hat[i] = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
  }
}

model {
  // Priors
  p_Am  ~ lognormal(prior_p_Am_mu, prior_p_Am_sd);
  p_M   ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa ~ beta(prior_kappa_a, prior_kappa_b);
  v     ~ lognormal(prior_v_mu, prior_v_sd);
  E_G   ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  E0    ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0    ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);

  // Likelihood
  L_obs ~ normal(L_hat, sigma_L);
}

generated quantities {
  array[N_obs] real L_rep;
  array[N_obs] real log_lik;

  for (i in 1:N_obs) {
    L_rep[i] = normal_rng(L_hat[i], sigma_L);
    log_lik[i] = normal_lpdf(L_obs[i] | L_hat[i], sigma_L);
  }
}

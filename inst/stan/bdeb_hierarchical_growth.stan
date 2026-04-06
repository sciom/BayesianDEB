// ============================================================
// BayesianDEB: Hierarchical multi-individual growth model
// Non-centred parameterisation for p_Am (individual variation)
// Shared: p_M, kappa, v, E_G across individuals
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
  int<lower=1> N_ind;
  int<lower=1> max_N_obs;
  array[N_ind] int<lower=1> N_obs;
  array[N_ind, max_N_obs] real t_obs;
  array[N_ind, max_N_obs] real L_obs;    // NaN for missing
  real<lower=0, upper=1> f_food;

  // Prior hyperparameters
  real prior_mu_log_p_Am_mu;
  real<lower=0> prior_mu_log_p_Am_sd;
  real<lower=0> prior_sigma_log_p_Am_rate;
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
  // Population-level
  real mu_log_p_Am;
  real<lower=0> sigma_log_p_Am;

  // Individual deviations (non-centred)
  array[N_ind] real z_log_p_Am;

  // Shared parameters
  real<lower=0> p_M;
  real<lower=0, upper=1> kappa;
  real<lower=0> v;
  real<lower=0> E_G;

  // Initial conditions (per individual)
  real<lower=0> E0;
  array[N_ind] real<lower=0> L0;

  // Observation error
  real<lower=0> sigma_L;
}

transformed parameters {
  array[N_ind] real<lower=0> p_Am_ind;

  for (j in 1:N_ind) {
    p_Am_ind[j] = exp(mu_log_p_Am + sigma_log_p_Am * z_log_p_Am[j]);
  }
}

model {
  // Hyperpriors
  mu_log_p_Am    ~ normal(prior_mu_log_p_Am_mu, prior_mu_log_p_Am_sd);
  sigma_log_p_Am ~ exponential(prior_sigma_log_p_Am_rate);

  // Individual deviations
  z_log_p_Am ~ std_normal();

  // Shared parameter priors
  p_M   ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa ~ beta(prior_kappa_a, prior_kappa_b);
  v     ~ lognormal(prior_v_mu, prior_v_sd);
  E_G   ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  E0    ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0    ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);

  // Per-individual ODE + likelihood
  for (j in 1:N_ind) {
    vector[2] x0;
    real V0 = L0[j] * L0[j] * L0[j];
    x0[1] = E0 * V0;
    x0[2] = V0;

    array[5] real theta = {p_Am_ind[j], p_M, kappa, v, E_G};
    array[1] real env = {f_food};

    // Collect valid observation times for this individual
    array[N_obs[j]] real t_ind;
    for (i in 1:N_obs[j]) {
      t_ind[i] = t_obs[j, i];
    }

    array[N_obs[j]] vector[2] x_sol = ode_bdf(
      deb_growth_ode, x0, 0.0, t_ind, theta, env, 1e-6, 1e-6, 1e4
    );

    for (i in 1:N_obs[j]) {
      real L_hat = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
      if (!is_nan(L_obs[j, i])) {
        L_obs[j, i] ~ normal(L_hat, sigma_L);
      }
    }
  }
}

generated quantities {
  // Population-level prediction for a new individual
  real p_Am_new = exp(normal_rng(mu_log_p_Am, sigma_log_p_Am));
}

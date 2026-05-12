// ============================================================
// BayesianDEB: Hierarchical multi-individual growth model
// Non-centred parameterisation for p_Am (individual variation)
// Growth obs: 1=normal, 2=lognormal, 3=student_t
// Within-chain parallelism via reduce_sum
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

  real partial_log_lik(array[] int ind_slice,
                       int start, int end,
                       int max_N_obs,
                       array[] int N_obs,
                       array[,] real t_obs,
                       array[,] real L_obs,
                       real f_food,
                       array[] real p_Am_ind,
                       real p_M, real kappa, real v, real E_G,
                       real E0, array[] real L0,
                       real sigma_L,
                       real c_T,
                       int obs_growth, real obs_nu) {
    real ll = 0;
    for (k in 1:size(ind_slice)) {
      int j = ind_slice[k];

      real V0 = L0[j] * L0[j] * L0[j];
      vector[2] x0 = [E0 * V0, V0]';

      array[N_obs[j]] real t_ind;
      for (i in 1:N_obs[j]) t_ind[i] = t_obs[j, i];

      array[N_obs[j]] vector[2] x_sol = ode_bdf_tol(
        deb_growth_ode, x0, 0.0, t_ind,
        1e-5, 1e-8, 1000000,
        p_Am_ind[j] * c_T, p_M * c_T, kappa, v * c_T, E_G, f_food
      );

      for (i in 1:N_obs[j]) {
        real L_hat = pow(fmax(x_sol[i][2], 1e-12), 1.0 / 3.0);
        if (!is_nan(L_obs[j, i])) {
          if (obs_growth == 1)
            ll += normal_lpdf(L_obs[j, i] | L_hat, sigma_L);
          else if (obs_growth == 2)
            ll += lognormal_lpdf(L_obs[j, i] | log(fmax(L_hat, 1e-12)), sigma_L);
          else
            ll += student_t_lpdf(L_obs[j, i] | obs_nu, L_hat, sigma_L);
        }
      }
    }
    return ll;
  }
}

data {
  int<lower=1> N_ind;
  int<lower=1> max_N_obs;
  array[N_ind] int<lower=1> N_obs;
  array[N_ind, max_N_obs] real t_obs;
  array[N_ind, max_N_obs] real L_obs;
  real<lower=0, upper=1> f_food;

  int<lower=1, upper=3> obs_growth;
  real<lower=1> obs_nu;
  int<lower=1, upper=2> obs_repro;   // unused here, declared for data compatibility

  int<lower=0, upper=1> has_temperature;
  real<lower=0> T_obs;
  real<lower=0> T_ref;
  real<lower=0> T_A;

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

transformed data {
  array[N_ind] int ind_idx;
  for (j in 1:N_ind) ind_idx[j] = j;

  real c_T = has_temperature
    ? exp(T_A / T_ref - T_A / T_obs)
    : 1.0;
}

parameters {
  real mu_log_p_Am;
  real<lower=0> sigma_log_p_Am;
  array[N_ind] real z_log_p_Am;
  real<lower=0> p_M;
  real<lower=0, upper=1> kappa;
  real<lower=0> v;
  real<lower=0> E_G;
  real<lower=0> E0;
  array[N_ind] real<lower=0> L0;
  real<lower=0> sigma_L;
}

transformed parameters {
  array[N_ind] real<lower=0> p_Am_ind;
  for (j in 1:N_ind) {
    p_Am_ind[j] = exp(mu_log_p_Am + sigma_log_p_Am * z_log_p_Am[j]);
  }
}

model {
  mu_log_p_Am    ~ normal(prior_mu_log_p_Am_mu, prior_mu_log_p_Am_sd);
  sigma_log_p_Am ~ exponential(prior_sigma_log_p_Am_rate);
  z_log_p_Am     ~ std_normal();
  p_M     ~ lognormal(prior_p_M_mu, prior_p_M_sd);
  kappa   ~ beta(prior_kappa_a, prior_kappa_b);
  v       ~ lognormal(prior_v_mu, prior_v_sd);
  E_G     ~ lognormal(prior_E_G_mu, prior_E_G_sd);
  E0      ~ lognormal(prior_E0_mu, prior_E0_sd);
  L0      ~ lognormal(prior_L0_mu, prior_L0_sd);
  sigma_L ~ normal(0, prior_sigma_L_sd);  // half-normal via <lower=0> constraint

  target += reduce_sum(partial_log_lik, ind_idx, 1,
                       max_N_obs, N_obs, t_obs, L_obs, f_food,
                       p_Am_ind, p_M, kappa, v, E_G,
                       E0, L0, sigma_L, c_T,
                       obs_growth, obs_nu);
}

generated quantities {
  real p_Am_new = exp(normal_rng(mu_log_p_Am, sigma_log_p_Am));
}

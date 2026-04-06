// ============================================================
// DEB ODE right-hand side functions
// Shared across all BayesianDEB Stan models
// ============================================================

// Standard DEB: 2-state growth model (reserve E, structure V)
// theta = [p_Am, p_M, kappa, v, E_G]
vector deb_growth_ode(real t, vector x, array[] real theta,
                      array[] real env) {
  real E = x[1];       // reserve (J)
  real V = x[2];       // structural volume (cm^3)

  real p_Am  = theta[1];
  real p_M   = theta[2];
  real kappa = theta[3];
  real v_    = theta[4];
  real E_G   = theta[5];
  real f     = env[1];  // scaled functional response

  real L = pow(V, 1.0 / 3.0);

  // Energy fluxes
  real p_A = f * p_Am * L * L;                     // assimilation
  real p_C = E * v_ * L / (E + E_G * V + 1e-12);  // mobilisation (per volume)
  real p_M_flux = p_M * V;                         // somatic maintenance

  // State derivatives
  real dE = p_A - p_C;
  real dV = (kappa * p_C - p_M_flux) / E_G;

  // Prevent structural shrinkage below zero
  if (V < 1e-12 && dV < 0) dV = 0;

  return [dE, dV]';
}

// Standard DEB: 3-state model (reserve E, structure V, repro buffer R)
// theta = [p_Am, p_M, kappa, v, E_G, k_J, E_Hp]
vector deb_growth_repro_ode(real t, vector x, array[] real theta,
                            array[] real env) {
  real E  = x[1];      // reserve (J)
  real V  = x[2];      // structural volume (cm^3)
  real R  = x[3];      // reproduction buffer (J)

  real p_Am  = theta[1];
  real p_M   = theta[2];
  real kappa = theta[3];
  real v_    = theta[4];
  real E_G   = theta[5];
  real k_J   = theta[6];
  real E_Hp  = theta[7];
  real f     = env[1];

  real L = pow(V, 1.0 / 3.0);

  // Energy fluxes
  real p_A = f * p_Am * L * L;
  real p_C = E * v_ * L / (E + E_G * V + 1e-12);
  real p_M_flux = p_M * V;
  real p_J = k_J * E_Hp;  // maturity maintenance (post-puberty)

  // Kappa-rule allocation
  real dE = p_A - p_C;
  real dV = (kappa * p_C - p_M_flux) / E_G;
  real dR = (1 - kappa) * p_C - p_J;

  if (V < 1e-12 && dV < 0) dV = 0;
  if (dR < 0) dR = 0;  // no negative reproduction flux

  return [dE, dV, dR]';
}

// DEBtox: 4-state model adding internal damage (scaled damage d)
// theta = [p_Am, p_M, kappa, v, E_G, k_d, z_w, b_w, mode]
// env   = [f, C_w]  where C_w = external concentration
// mode: 1=assimilation stress, 2=maintenance stress, 3=growth cost stress
vector deb_tox_ode(real t, vector x, array[] real theta,
                   array[] real env) {
  real E  = x[1];
  real V  = x[2];
  real R  = x[3];
  real Dw = x[4];      // scaled internal damage

  real p_Am  = theta[1];
  real p_M   = theta[2];
  real kappa = theta[3];
  real v_    = theta[4];
  real E_G   = theta[5];
  real k_d   = theta[6];  // damage recovery rate
  real z_w   = theta[7];  // no-effect concentration (NEC)
  real b_w   = theta[8];  // killing rate / effect intensity
  int  mode  = 1;         // mode of action (cast from real below)

  real f   = env[1];
  real C_w = env[2];      // external concentration

  // Toxicokinetics: scaled damage
  real dDw = k_d * (fmax(C_w - z_w, 0.0) - Dw);

  // Stress factor
  real s = b_w * fmax(Dw, 0.0);

  real L = pow(V, 1.0 / 3.0);

  // Base fluxes with stress
  real p_A;
  real p_M_flux;
  real eff_E_G;

  // Mode of action: stress on assimilation
  p_A      = f * p_Am * L * L * fmax(1.0 - s, 0.0);
  p_M_flux = p_M * V;
  eff_E_G  = E_G;

  real p_C = E * v_ * L / (E + eff_E_G * V + 1e-12);

  real dE = p_A - p_C;
  real dV = (kappa * p_C - p_M_flux) / eff_E_G;
  real dR = fmax((1 - kappa) * p_C - 0.0, 0.0);

  if (V < 1e-12 && dV < 0) dV = 0;

  return [dE, dV, dR, dDw]';
}

// Arrhenius temperature correction factor
real arrhenius_correction(real T, real T_ref, real T_A) {
  return exp(T_A / T_ref - T_A / T);
}

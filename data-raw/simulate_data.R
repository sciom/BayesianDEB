# =============================================================
# Simulate example datasets for BayesianDEB package
# Based on BDEB book case studies (Ch. 20-21)
# =============================================================

set.seed(42)

# --- Helper: simple DEB growth ODE (Euler integration) ---
sim_deb_growth <- function(t_max, dt, pars, f = 1.0) {
	n <- ceiling(t_max / dt)
	t <- seq(0, t_max, length.out = n + 1)
	E <- numeric(n + 1)
	V <- numeric(n + 1)

	V[1] <- pars$L0^3
	E[1] <- pars$E0 * V[1]

	for (i in seq_len(n)) {
		L <- V[i]^(1/3)
		p_A <- f * pars$p_Am * L^2
		p_C <- E[i] * pars$v * L / (E[i] + pars$E_G * V[i] + 1e-12)
		p_M <- pars$p_M * V[i]

		dE <- p_A - p_C
		dV <- (pars$kappa * p_C - p_M) / pars$E_G

		E[i + 1] <- max(E[i] + dE * dt, 1e-6)
		V[i + 1] <- max(V[i] + dV * dt, 1e-6)
	}
	data.frame(time = t, E = E, V = V, L = V^(1/3))
}

# ============================================================
# 1. eisenia_growth: Eisenia fetida growth data
#    21 individuals, 12 weeks, weekly measurements
# ============================================================
true_pars_eisenia <- list(
	p_Am  = 5.0,    # J/d/cm^2
	p_M   = 0.5,    # J/d/cm^3
	kappa = 0.75,
	v     = 0.2,    # cm/d
	E_G   = 400,    # J/cm^3
	E0    = 1.0,    # J/cm^3
	L0    = 0.1     # cm
)

sigma_L <- 0.015
t_obs_weeks <- 0:12
t_obs <- t_obs_weeks * 7  # convert to days

n_ind <- 21
eisenia_growth_list <- list()

for (j in seq_len(n_ind)) {
	# Small individual variation in p_Am
	ind_pars <- true_pars_eisenia
	ind_pars$p_Am <- rlnorm(1, log(true_pars_eisenia$p_Am), 0.1)
	ind_pars$L0   <- rlnorm(1, log(true_pars_eisenia$L0), 0.05)

	traj <- sim_deb_growth(max(t_obs), dt = 0.1, ind_pars, f = 1.0)

	# Interpolate to observation times
	L_true <- approx(traj$time, traj$L, xout = t_obs)$y
	L_obs  <- L_true + rnorm(length(t_obs), 0, sigma_L)
	L_obs  <- pmax(L_obs, 0.01)

	eisenia_growth_list[[j]] <- data.frame(
		id     = j,
		time   = t_obs,
		length = round(L_obs, 4)
	)
}

eisenia_growth <- do.call(rbind, eisenia_growth_list)

# ============================================================
# 2. folsomia_repro: Folsomia candida reproduction test
#    5 Cd concentrations, 6 replicates each, 28 days
# ============================================================
concentrations <- c(0, 10, 50, 100, 500)
n_reps <- 6
t_harvest <- 28

folsomia_repro_list <- list()
id_counter <- 0

for (conc in concentrations) {
	# Effect: reduce reproduction with increasing concentration
	nec <- 30  # no-effect concentration
	b_w <- 0.005
	stress <- b_w * max(conc - nec, 0)
	base_offspring <- 120
	expected <- base_offspring * max(1 - stress, 0)

	for (rep in seq_len(n_reps)) {
		id_counter <- id_counter + 1
		count <- rpois(1, lambda = max(expected + rnorm(1, 0, 10), 0))

		folsomia_repro_list[[id_counter]] <- data.frame(
			id            = id_counter,
			concentration = conc,
			t_start       = 0,
			t_end         = t_harvest,
			count         = count
		)
	}
}

folsomia_repro <- do.call(rbind, folsomia_repro_list)

# ============================================================
# 3. debtox_growth: Growth under toxicant exposure
#    4 concentration groups, 10 individuals each, 6 weeks
# ============================================================
conc_groups <- c(0, 20, 80, 200)
n_per_group <- 10
t_obs_tox <- seq(0, 42, by = 7)

debtox_growth_list <- list()
id_counter <- 0

for (cg in conc_groups) {
	nec_true <- 15
	b_w_true <- 0.003
	stress <- b_w_true * max(cg - nec_true, 0)

	for (rep in seq_len(n_per_group)) {
		id_counter <- id_counter + 1
		ind_pars <- true_pars_eisenia
		# Stress reduces effective assimilation
		ind_pars$p_Am <- true_pars_eisenia$p_Am * max(1 - stress, 0)
		ind_pars$p_Am <- ind_pars$p_Am * rlnorm(1, 0, 0.05)

		traj <- sim_deb_growth(max(t_obs_tox), dt = 0.1, ind_pars, f = 1.0)
		L_true <- approx(traj$time, traj$L, xout = t_obs_tox)$y
		L_obs  <- L_true + rnorm(length(t_obs_tox), 0, sigma_L)
		L_obs  <- pmax(L_obs, 0.01)

		debtox_growth_list[[id_counter]] <- data.frame(
			id            = id_counter,
			concentration = cg,
			time          = t_obs_tox,
			length        = round(L_obs, 4)
		)
	}
}

debtox_growth <- do.call(rbind, debtox_growth_list)

# ============================================================
# Save datasets
# ============================================================
usethis::use_data(eisenia_growth, overwrite = TRUE)
usethis::use_data(folsomia_repro, overwrite = TRUE)
usethis::use_data(debtox_growth, overwrite = TRUE)

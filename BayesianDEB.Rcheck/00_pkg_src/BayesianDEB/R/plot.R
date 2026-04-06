#' Plot Methods for BDEB Objects
#'
#' Visualisation methods for BDEB fits, posterior predictive checks,
#' and derived quantities using ggplot2.  All plots return ggplot2
#' objects that can be further customised.
#'
#' @importFrom rlang .data
#' @name bdeb_plots
NULL

#' Plot a BDEB Fit
#'
#' @param x A [bdeb_fit()] object.
#' @param type Type of plot. One of:
#'   - `"trace"`: MCMC trace plots
#'   - `"posterior"`: marginal posterior densities
#'   - `"pairs"`: bivariate posterior scatter plots
#'   - `"trajectory"`: predicted trajectories with data overlay
#' @param pars Character vector of parameters to plot. Default: core DEB
#'   parameters.
#' @param n_draws Number of posterior draws for trajectory plots. Default 100.
#' @param seed Integer seed for reproducible draw selection in trajectory
#'   plots.  Default `NULL` (no seed).
#' @param ... Additional arguments passed to bayesplot functions.
#' @return A ggplot2 object.
#' @export
plot.bdeb_fit <- function(x, type = c("trace", "posterior", "pairs",
                                       "trajectory"),
                          pars = NULL, n_draws = 100, seed = NULL, ...) {
	type <- match.arg(type)

	if (is.null(pars)) {
		pars <- get_core_pars(x$model$type)
	}

	draws <- x$fit$draws()

	switch(type,
		trace     = plot_trace(draws, pars, ...),
		posterior = plot_posterior(draws, pars, ...),
		pairs     = plot_pairs(draws, pars, ...),
		trajectory = plot_trajectory(x, n_draws, seed, ...)
	)
}

#' Plot Posterior Predictive Checks
#'
#' @param x A [bdeb_ppc()] object.
#' @param n_draws Number of replicated trajectories to show. Default 50.
#' @param ... Ignored.
#' @return A ggplot2 object.
#' @export
plot.bdeb_ppc <- function(x, n_draws = 50, ...) {
	if (!is.null(x$growth)) {
		plot_ppc_growth(x$growth, n_draws)
	} else if (!is.null(x$reproduction)) {
		plot_ppc_repro(x$reproduction, n_draws)
	} else {
		cli::cli_abort("No PPC data available to plot.")
	}
}

# --- Internal plotting functions ---

#' @keywords internal
get_core_pars <- function(model_type) {
	switch(model_type,
		individual   = c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"),
		growth_repro = c("p_Am", "p_M", "kappa", "v", "E_G", "k_J", "sigma_L", "k_R"),
		hierarchical = c("mu_log_p_Am", "sigma_log_p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"),
		debtox       = c("p_Am", "p_M", "kappa", "k_d", "z_w", "b_w", "sigma_L"),
		c("p_Am", "p_M", "kappa", "sigma_L")
	)
}

#' @keywords internal
plot_trace <- function(draws, pars, ...) {
	bayesplot::mcmc_trace(draws, pars = pars, ...) +
		ggplot2::theme_bw() +
		ggplot2::labs(title = "MCMC Trace Plots")
}

#' @keywords internal
plot_posterior <- function(draws, pars, ...) {
	bayesplot::mcmc_dens_overlay(draws, pars = pars, ...) +
		ggplot2::theme_bw() +
		ggplot2::labs(title = "Posterior Densities")
}

#' @keywords internal
plot_pairs <- function(draws, pars, ...) {
	bayesplot::mcmc_pairs(draws, pars = pars, ...) +
		ggplot2::labs(title = "Posterior Pairs")
}

# --- Trajectory: dispatcher ---

#' @keywords internal
plot_trajectory <- function(fit, n_draws, seed = NULL, ...) {
	if (!is.null(seed)) set.seed(seed)
	switch(fit$model$type,
		individual   = plot_trajectory_individual(fit, n_draws),
		growth_repro = plot_trajectory_individual(fit, n_draws),
		hierarchical = plot_trajectory_hierarchical(fit, n_draws),
		debtox       = plot_trajectory_debtox(fit, n_draws)
	)
}

# --- Trajectory: individual / growth_repro ---
# These models store L_hat in transformed parameters

#' @keywords internal
plot_trajectory_individual <- function(fit, n_draws) {
	draws <- posterior::as_draws_df(fit$fit$draws())
	L_hat_vars <- grep("^L_hat\\[", names(draws), value = TRUE)

	if (length(L_hat_vars) == 0) {
		cli::cli_abort("No {.var L_hat} variables found in posterior draws.")
	}

	L_hat <- as.matrix(draws[, L_hat_vars])
	n_total <- nrow(L_hat)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	t_obs <- fit$model$stan_data$t_obs
	L_obs <- fit$model$stan_data$L_obs
	if (fit$model$type == "growth_repro") {
		t_obs <- fit$model$stan_data$t_L
	}

	traj_list <- lapply(idx, function(i) {
		data.frame(
			time   = t_obs[seq_len(ncol(L_hat))],
			length = as.numeric(L_hat[i, ]),
			draw   = i
		)
	})
	traj_df <- do.call(rbind, traj_list)
	obs_df <- data.frame(time = t_obs[seq_along(L_obs)], length = L_obs)

	ggplot2::ggplot() +
		ggplot2::geom_line(
			data = traj_df,
			ggplot2::aes(x = .data$time, y = .data$length, group = .data$draw),
			alpha = 0.15, colour = "steelblue"
		) +
		ggplot2::geom_point(
			data = obs_df,
			ggplot2::aes(x = .data$time, y = .data$length),
			size = 2.5, colour = "black"
		) +
		ggplot2::theme_bw() +
		ggplot2::labs(title = "Posterior Predicted Trajectories",
		             x = "Time", y = "Structural Length")
}

# --- Trajectory: hierarchical ---
# ODE is solved inside reduce_sum so no L_hat is stored.
# We forward-simulate from posterior draws in R (Euler).

#' @keywords internal
plot_trajectory_hierarchical <- function(fit, n_draws) {
	draws <- posterior::as_draws_df(fit$fit$draws())
	sd <- fit$model$stan_data
	n_ind <- sd$N_ind
	n_total <- nrow(draws)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	# Collect all observed data for all individuals
	obs_list <- list()
	for (j in seq_len(n_ind)) {
		ni <- sd$N_obs[j]
		tt <- sd$t_obs[j, seq_len(ni)]
		ll <- sd$L_obs[j, seq_len(ni)]
		valid <- !is.nan(ll)
		obs_list[[j]] <- data.frame(
			time = tt[valid], length = ll[valid],
			individual = factor(j)
		)
	}
	obs_df <- do.call(rbind, obs_list)

	# Simulate trajectories for each draw x each individual
	traj_list <- list()
	for (i in idx) {
		p_M   <- draws$p_M[i]
		kappa <- draws$kappa[i]
		v_val <- draws$v[i]
		E_G   <- draws$E_G[i]
		E0    <- draws$E0[i]
		sigma <- draws$sigma_L[i]
		f     <- sd$f_food

		for (j in seq_len(n_ind)) {
			p_Am_j <- draws[[paste0("p_Am_ind[", j, "]")]][i]
			L0_j   <- draws[[paste0("L0[", j, "]")]][i]
			ni     <- sd$N_obs[j]
			t_end  <- sd$t_obs[j, ni]

			traj <- sim_deb_euler(t_end, p_Am_j, p_M, kappa, v_val,
			                      E_G, E0, L0_j, f)
			t_out <- sd$t_obs[j, seq_len(ni)]
			L_pred <- stats::approx(traj$time, traj$L, xout = t_out)$y

			traj_list[[length(traj_list) + 1L]] <- data.frame(
				time = t_out, length = L_pred,
				individual = factor(j), draw = i
			)
		}
	}
	traj_df <- do.call(rbind, traj_list)

	ggplot2::ggplot() +
		ggplot2::geom_line(
			data = traj_df,
			ggplot2::aes(x = .data$time, y = .data$length,
			             group = interaction(.data$draw, .data$individual)),
			alpha = 0.08, colour = "steelblue"
		) +
		ggplot2::geom_point(
			data = obs_df,
			ggplot2::aes(x = .data$time, y = .data$length),
			size = 1.5, colour = "black"
		) +
		ggplot2::facet_wrap(~individual, scales = "free_y") +
		ggplot2::theme_bw(base_size = 10) +
		ggplot2::labs(
			title = "Hierarchical Model: Posterior Trajectories by Individual",
			x = "Time", y = "Structural Length"
		)
}

# --- Trajectory: debtox ---
# ODE is solved per group inside reduce_sum.
# We forward-simulate with toxicant stress in R.

#' @keywords internal
plot_trajectory_debtox <- function(fit, n_draws) {
	draws <- posterior::as_draws_df(fit$fit$draws())
	sd <- fit$model$stan_data
	n_groups <- sd$N_groups
	n_total <- nrow(draws)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	# Observed data per group
	obs_list <- list()
	for (g in seq_len(n_groups)) {
		ni <- sd$N_obs[g]
		tt <- sd$t_obs[g, seq_len(ni)]
		ll <- sd$L_obs[g, seq_len(ni)]
		valid <- !is.nan(ll)
		obs_list[[g]] <- data.frame(
			time = tt[valid], length = ll[valid],
			group = factor(paste0("C = ", sd$C_w[g]))
		)
	}
	obs_df <- do.call(rbind, obs_list)

	# Simulate trajectories per draw x group
	traj_list <- list()
	for (i in idx) {
		p_Am  <- draws$p_Am[i]
		p_M   <- draws$p_M[i]
		kappa <- draws$kappa[i]
		v_val <- draws$v[i]
		E_G   <- draws$E_G[i]
		E0    <- draws$E0[i]
		L0    <- draws$L0[i]
		k_d   <- draws$k_d[i]
		z_w   <- draws$z_w[i]
		b_w   <- draws$b_w[i]
		f     <- sd$f_food

		for (g in seq_len(n_groups)) {
			C_w <- sd$C_w[g]
			ni  <- sd$N_obs[g]
			t_end <- sd$t_obs[g, ni]

			traj <- sim_debtox_euler(t_end, p_Am, p_M, kappa, v_val,
			                         E_G, E0, L0, f, k_d, z_w, b_w, C_w)
			t_out <- sd$t_obs[g, seq_len(ni)]
			L_pred <- stats::approx(traj$time, traj$L, xout = t_out)$y

			traj_list[[length(traj_list) + 1L]] <- data.frame(
				time = t_out, length = L_pred,
				group = factor(paste0("C = ", C_w)), draw = i
			)
		}
	}
	traj_df <- do.call(rbind, traj_list)

	ggplot2::ggplot() +
		ggplot2::geom_line(
			data = traj_df,
			ggplot2::aes(x = .data$time, y = .data$length,
			             group = .data$draw),
			alpha = 0.1, colour = "steelblue"
		) +
		ggplot2::geom_point(
			data = obs_df,
			ggplot2::aes(x = .data$time, y = .data$length),
			size = 2, colour = "black"
		) +
		ggplot2::facet_wrap(~group) +
		ggplot2::theme_bw(base_size = 11) +
		ggplot2::labs(
			title = "DEBtox: Posterior Trajectories by Concentration Group",
			x = "Time", y = "Structural Length"
		)
}

# --- Euler forward simulators for R-side trajectory plots ---

#' @keywords internal
sim_deb_euler <- function(t_max, p_Am, p_M, kappa, v, E_G, E0, L0, f,
                          dt = 0.5) {
	n <- ceiling(t_max / dt)
	t <- seq(0, t_max, length.out = n + 1)
	E <- V <- numeric(n + 1)
	V[1] <- L0^3
	E[1] <- E0 * V[1]

	for (i in seq_len(n)) {
		L <- V[i]^(1/3)
		pA <- f * p_Am * L^2
		pC <- E[i] * v * L / (E[i] + E_G * V[i] + 1e-12)  # consistent with Stan ODE
		pM <- p_M * V[i]
		dE <- pA - pC
		dV <- (kappa * pC - pM) / E_G
		E[i + 1] <- max(E[i] + dE * dt, 1e-12)
		V[i + 1] <- max(V[i] + dV * dt, 1e-12)
	}
	data.frame(time = t, L = V^(1/3))
}

#' @keywords internal
sim_debtox_euler <- function(t_max, p_Am, p_M, kappa, v, E_G, E0, L0,
                             f, k_d, z_w, b_w, C_w, dt = 0.5) {
	n <- ceiling(t_max / dt)
	t <- seq(0, t_max, length.out = n + 1)
	E <- V <- Dw <- numeric(n + 1)
	V[1] <- L0^3
	E[1] <- E0 * V[1]
	Dw[1] <- 0

	for (i in seq_len(n)) {
		L <- V[i]^(1/3)
		dDw <- k_d * (max(C_w - z_w, 0) - Dw[i])
		s <- b_w * max(Dw[i], 0)
		pA <- f * p_Am * L^2 * max(1 - s, 0)
		pC <- E[i] * v * L / (E[i] + E_G * V[i] + 1e-12)  # consistent with Stan ODE
		pM <- p_M * V[i]
		dE <- pA - pC
		dV <- (kappa * pC - pM) / E_G
		E[i + 1]  <- max(E[i] + dE * dt, 1e-12)
		V[i + 1]  <- max(V[i] + dV * dt, 1e-12)
		Dw[i + 1] <- max(Dw[i] + dDw * dt, 0)
	}
	data.frame(time = t, L = V^(1/3))
}

# --- PPC plotting ---

#' @keywords internal
plot_ppc_growth <- function(growth, n_draws) {
	L_rep <- growth$L_rep
	L_obs <- growth$L_obs
	t_obs <- growth$t_obs

	n_total <- nrow(L_rep)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	traj_list <- lapply(idx, function(i) {
		data.frame(time = t_obs, length = as.numeric(L_rep[i, ]), draw = i)
	})
	traj_df <- do.call(rbind, traj_list)
	obs_df <- data.frame(time = t_obs, length = L_obs)

	ggplot2::ggplot() +
		ggplot2::geom_line(
			data = traj_df,
			ggplot2::aes(x = .data$time, y = .data$length, group = .data$draw),
			alpha = 0.15, colour = "grey60"
		) +
		ggplot2::geom_point(
			data = obs_df,
			ggplot2::aes(x = .data$time, y = .data$length),
			size = 2.5, colour = "red"
		) +
		ggplot2::theme_bw() +
		ggplot2::labs(title = "Posterior Predictive Check: Growth",
		             x = "Time", y = "Structural Length")
}

#' @keywords internal
plot_ppc_repro <- function(repro, n_draws) {
	R_rep <- repro$R_rep
	R_obs <- repro$R_obs

	n_total <- nrow(R_rep)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	rep_df <- data.frame(count = as.numeric(R_rep[idx, ]), type = "Predicted")
	obs_df <- data.frame(count = R_obs, type = "Observed")

	ggplot2::ggplot() +
		ggplot2::geom_histogram(
			data = rep_df,
			ggplot2::aes(x = .data$count),
			fill = "steelblue", alpha = 0.4, bins = 30
		) +
		ggplot2::geom_vline(
			xintercept = R_obs,
			colour = "red", linewidth = 0.8
		) +
		ggplot2::theme_bw() +
		ggplot2::labs(title = "Posterior Predictive Check: Reproduction",
		             x = "Offspring Count", y = "Frequency")
}

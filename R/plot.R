#' Plot Methods for BDEB Objects
#'
#' Publication-quality visualisation of BDEB fits, posterior predictive
#' checks, and derived quantities using ggplot2.
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
#' @param ... Additional arguments passed to bayesplot functions.
#' @return A ggplot2 object.
#' @export
plot.bdeb_fit <- function(x, type = c("trace", "posterior", "pairs",
                                       "trajectory"),
                          pars = NULL, n_draws = 100, ...) {
	type <- match.arg(type)

	if (is.null(pars)) {
		pars <- get_core_pars(x$model$type)
	}

	draws <- x$fit$draws()

	switch(type,
		trace     = plot_trace(draws, pars, ...),
		posterior = plot_posterior(draws, pars, ...),
		pairs     = plot_pairs(draws, pars, ...),
		trajectory = plot_trajectory(x, n_draws, ...)
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

#' @keywords internal
plot_trajectory <- function(fit, n_draws, ...) {
	draws <- posterior::as_draws_df(fit$fit$draws())
	L_hat_vars <- grep("^L_hat\\[", names(draws), value = TRUE)

	if (length(L_hat_vars) == 0) {
		cli::cli_abort("No L_hat variables found for trajectory plot.")
	}

	L_hat <- as.matrix(draws[, L_hat_vars])
	n_total <- nrow(L_hat)
	idx <- seq(1, n_total, length.out = min(n_draws, n_total))
	idx <- round(idx)

	# Get times
	t_obs <- fit$model$stan_data$t_obs
	if (is.matrix(t_obs)) t_obs <- t_obs[1, ]
	n_t <- length(t_obs)

	# Build data frame for ggplot
	traj_list <- lapply(idx, function(i) {
		data.frame(
			time   = t_obs[seq_len(ncol(L_hat))],
			length = as.numeric(L_hat[i, ]),
			draw   = i
		)
	})
	traj_df <- do.call(rbind, traj_list)

	# Observed data
	L_obs <- fit$model$stan_data$L_obs
	if (is.matrix(L_obs)) L_obs <- L_obs[1, ]
	obs_df <- data.frame(
		time   = t_obs[seq_along(L_obs)],
		length = L_obs
	)

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
		ggplot2::labs(
			title = "Posterior Predicted Trajectories",
			x = "Time",
			y = "Structural Length"
		)
}

#' @keywords internal
plot_ppc_growth <- function(growth, n_draws) {
	L_rep <- growth$L_rep
	L_obs <- growth$L_obs
	t_obs <- growth$t_obs

	n_total <- nrow(L_rep)
	idx <- seq(1, n_total, length.out = min(n_draws, n_total))
	idx <- round(idx)

	traj_list <- lapply(idx, function(i) {
		data.frame(
			time   = t_obs,
			length = as.numeric(L_rep[i, ]),
			draw   = i
		)
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
		ggplot2::labs(
			title = "Posterior Predictive Check: Growth",
			x = "Time",
			y = "Structural Length"
		)
}

#' @keywords internal
plot_ppc_repro <- function(repro, n_draws) {
	R_rep <- repro$R_rep
	R_obs <- repro$R_obs

	n_total <- nrow(R_rep)
	idx <- seq(1, n_total, length.out = min(n_draws, n_total))
	idx <- round(idx)

	# Simple histogram overlay
	rep_df <- data.frame(
		count = as.numeric(R_rep[idx, ]),
		type  = "Predicted"
	)
	obs_df <- data.frame(
		count = R_obs,
		type  = "Observed"
	)

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
		ggplot2::labs(
			title = "Posterior Predictive Check: Reproduction",
			x = "Offspring Count",
			y = "Frequency"
		)
}

#' DEBtox Model Specification
#'
#' Convenience wrapper for [bdeb_model()] that sets `type = "debtox"` and
#' provides a `stress` argument to select the physiological mode of action
#' (PMoA) of the toxicant.  The underlying TKTD framework follows
#' Jager et al. (2006) and the GUTS-RED-SD simplification of
#' Jager & Zimmer (2012):
#'
#' **Toxicokinetics.** Scaled internal damage \eqn{D_w} tracks the
#' external concentration with first-order kinetics:
#' \deqn{\frac{dD_w}{dt} = k_d\bigl(\max(C_w - z_w,\, 0) - D_w\bigr)}
#' where \eqn{k_d} is the dominant rate constant, \eqn{z_w} is the NEC
#' (no-effect concentration), and \eqn{C_w} is the external concentration.
#'
#' **Stress on assimilation.** The assimilation flux is reduced by a
#' factor \eqn{\max(1 - b_w D_w, 0)}, where \eqn{b_w} is the effect
#' intensity.  At steady state (\eqn{D_w = C_w - z_w}), the
#' \eqn{\mathrm{EC}_{50}} for 50\% assimilation reduction is
#' \eqn{z_w + 0.5/b_w}.
#'
#' @param data A [bdeb_data()] object with `concentration` specified.
#' @param stress Physiological mode of action.  Currently only
#'   `"assimilation"` is implemented.
#' @param priors Named list of priors.  Missing entries filled from
#'   `prior_default("debtox")`.  The toxicological parameters \eqn{k_d},
#'   \eqn{z_w}, \eqn{b_w} default to weakly informative log-normal priors.
#' @param ... Additional arguments passed to [bdeb_model()].
#' @return A `bdeb_model` object of type `"debtox"`.
#'
#' @note **Lifecycle: beta.**  Stress on assimilation is implemented and
#'   tested; `"maintenance"` and `"growth_cost"` modes are planned.
#'
#'   **Important limitation:** the DEBtox model fits **one ODE per
#'   concentration group**, not per individual.  If your data contain
#'   multiple individuals per concentration, they are automatically
#'   aggregated to group means per time point (with a warning).
#'   This is appropriate for group-level summary data but does not
#'   capture within-group individual variation.  A hierarchical
#'   DEBtox extension is planned for a future version.
#'
#' @references
#' Jager, T., Heugens, E.H.W. and Kooijman, S.A.L.M. (2006). Making
#' sense of ecotoxicological test results: towards application of
#' process-based models. *Ecotoxicology*, 15(3), 305--314.
#' \doi{10.1007/s10646-006-0060-x}
#'
#' Jager, T. and Zimmer, E.I. (2012). Simplified Dynamic Energy Budget
#' model for analysing ecotoxicity data. *Ecological Modelling*, 225,
#' 74--81. \doi{10.1016/j.ecolmodel.2011.11.012}
#' @export
#' @examples
#' # R-side specification only (no Stan sampling)
#' data(debtox_growth)
#' # one replicate per concentration avoids aggregation warning
#' dt <- debtox_growth[debtox_growth$id %in% c(1, 11, 21, 31), ]
#' dat <- bdeb_data(growth = dt, concentration = c(0, 20, 80, 200))
#' mod <- bdeb_tox(dat, stress = "assimilation")
bdeb_tox <- function(data,
                     stress = c("assimilation", "maintenance", "growth_cost"),
                     priors = list(),
                     ...) {
	stress <- match.arg(stress)

	if (is.null(data$concentration)) {
		cli::cli_abort("DEBtox requires {.arg concentration} in {.fun bdeb_data}.")
	}

	if (stress != "assimilation") {
		cli::cli_alert_warning(
			"Only 'assimilation' stress mode is currently implemented. Using assimilation."
		)
		stress <- "assimilation"
	}

	bdeb_model(data, type = "debtox", priors = priors, ...)
}

#' Extract EC50 and NEC from a DEBtox Fit
#'
#' Extracts the full posterior distribution of the \eqn{\mathrm{EC}_{50}}
#' and the NEC (no-effect concentration, \eqn{z_w}) from a fitted DEBtox
#' model.  Both quantities are computed analytically in the Stan
#' `generated quantities` block, avoiding the need for post-hoc root
#' finding.  At toxicokinetic steady state the stress factor equals
#' \eqn{s = b_w(C_w - z_w)} for \eqn{C_w > z_w}, so setting \eqn{s = 0.5}
#' yields
#'
#' \deqn{\mathrm{EC}_{50} = z_w + \frac{0.5}{b_w}.}
#'
#' The NEC is the threshold concentration below which no effect occurs;
#' it corresponds directly to the parameter \eqn{z_w} in the damage
#' model of Kooijman & Bedaux (1996).
#'
#' @references
#' Kooijman, S.A.L.M. and Bedaux, J.J.M. (1996). *The Analysis of
#' Aquatic Toxicity Data*. VU University Press, Amsterdam.
#'
#' @param fit A [bdeb_fit()] object from a DEBtox model.
#' @param prob Credible interval probability. Default 0.90.
#' @param verbose Logical; if `TRUE` (default) the summary table is
#'   printed via [cli::cli_verbatim()] / [message()] and can be
#'   silenced with [suppressMessages()].  Set to `FALSE` for a silent
#'   run; the invisible return value is identical.
#' @return A named list with:
#'   - `draws`: posterior draws of EC50
#'   - `summary`: mean, median, sd, lower, upper
#'   - `NEC`: posterior summary of the no-effect concentration
#' @export
bdeb_ec50 <- function(fit, prob = 0.90, verbose = TRUE) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	if (fit$model$type != "debtox") {
		cli::cli_abort("EC50 extraction requires a DEBtox model fit.")
	}

	draws <- posterior::as_draws_df(fit$fit$draws())
	alpha <- (1 - prob) / 2

	ec50_draws <- draws$EC50
	nec_draws  <- draws$NEC

	if (is.null(ec50_draws) || is.null(nec_draws)) {
		cli::cli_abort(c(
			"EC50 and/or NEC variables not found in posterior draws.",
			"i" = "Ensure the model was fitted with {.fn bdeb_tox}."
		))
	}

	ec50_summary <- data.frame(
		parameter = "EC50",
		mean   = mean(ec50_draws),
		median = stats::median(ec50_draws),
		sd     = stats::sd(ec50_draws),
		lower  = stats::quantile(ec50_draws, alpha),
		upper  = stats::quantile(ec50_draws, 1 - alpha)
	)

	nec_summary <- data.frame(
		parameter = "NEC",
		mean   = mean(nec_draws),
		median = stats::median(nec_draws),
		sd     = stats::sd(nec_draws),
		lower  = stats::quantile(nec_draws, alpha),
		upper  = stats::quantile(nec_draws, 1 - alpha)
	)

	result <- list(
		draws   = ec50_draws,
		summary = rbind(ec50_summary, nec_summary),
		NEC     = as.numeric(nec_draws)
	)

	if (verbose) {
		cli::cli_h3("DEBtox Effect Concentrations")
		# Route table through message() (CRAN-suppressible).
		tbl_lines <- utils::capture.output(
			print(result$summary, row.names = FALSE, digits = 3)
		)
		cli::cli_verbatim(tbl_lines)
	}

	invisible(result)
}

#' Plot DEBtox Dose-Response
#'
#' Produces a dose-response curve by **forward-simulating** the full
#' 4-state DEBtox ODE from the posterior in R (Euler integration).
#' For each posterior draw, the ODE is solved at every concentration in
#' a fine grid from 0 to \eqn{1.2 \times \max(C_w)}.  The y-axis shows
#' the predicted final structural length relative to the control (C = 0)
#' at the same draw, so that each curve propagates the full parameter
#' uncertainty through the dynamic model.
#'
#' This is a **visualisation tool**, not exact Stan inference.  The
#' R-side Euler integrator (step size `dt`) is an approximation of the
#' BDF solver used during fitting.  For quantitative
#' results, use [bdeb_ec50()] which extracts the analytically computed
#' EC\eqn{_{50}} and NEC directly from the Stan posterior.
#'
#' **Performance note.** Each draw requires `n_conc` ODE integrations,
#' so `n_draws * n_conc` total.  With default settings (100 draws,
#' 50 concentrations) this takes a few seconds.  Reduce `n_draws` for
#' faster interactive use.
#'
#' @param fit A [bdeb_fit()] object from a DEBtox model.
#' @param endpoint Which endpoint to plot.  Currently only `"growth"`.
#' @param n_draws Number of posterior draws to use.  Default 100.
#' @param n_conc Number of concentration points in the continuous grid.
#'   Default 50.
#' @param dt Euler integration step size (days).  Default 1.0.  Smaller
#'   values are more accurate but slower.  The Stan model uses BDF with
#'   adaptive stepping; this is an approximation for visualisation.
#' @param t_end End time for the simulation (days).  Default `NULL`
#'   (uses the last observation time from the fitted data).
#' @param seed Integer seed for reproducible draw selection.
#'   Default `NULL`.
#' @return A ggplot2 object.
#' @export
plot_dose_response <- function(fit, endpoint = "growth", n_draws = 100,
                               n_conc = 50, dt = 1.0, t_end = NULL,
                               seed = NULL) {
	if (!inherits(fit, "bdeb_fit") || fit$model$type != "debtox") {
		cli::cli_abort("Requires a fitted DEBtox model.")
	}

	draws <- posterior::as_draws_df(fit$fit$draws())
	n_total <- nrow(draws)
	if (!is.null(seed)) set.seed(seed)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	sd <- fit$model$stan_data
	C_w <- sd$C_w
	N_groups <- length(C_w)
	N_obs <- sd$N_obs

	# Final observation time
	if (is.null(t_end)) {
		t_end <- max(vapply(seq_len(N_groups), function(g) {
			sd$t_obs[g, N_obs[g]]
		}, numeric(1)))
	}

	# Concentration grid for the continuous curve
	c_seq <- seq(0, max(C_w) * 1.2, length.out = n_conc)

	# Observed: mean final length per group, normalised to control
	obs_final <- vapply(seq_len(N_groups), function(g) {
		vals <- sd$L_obs[g, seq_len(N_obs[g])]
		vals <- vals[!is.nan(vals)]
		if (length(vals) > 0) utils::tail(vals, 1) else NA_real_
	}, numeric(1))
	ctrl_idx <- which.min(C_w)
	ctrl_obs <- obs_final[ctrl_idx]
	if (min(C_w) > 0) {
		cli::cli_warn("No zero-concentration control group found; normalising to lowest concentration ({C_w[ctrl_idx]}).")
	}
	if (!is.finite(ctrl_obs) || ctrl_obs < 1e-12) ctrl_obs <- NA_real_
	obs_df <- data.frame(
		concentration = C_w,
		relative      = obs_final / ctrl_obs
	)

	# For each posterior draw, simulate the full ODE at each concentration
	pred_list <- lapply(idx, function(i) {
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

		# Final length at each concentration via full ODE
		L_final <- vapply(c_seq, function(cw) {
			traj <- sim_debtox_lsoda(t_end, p_Am, p_M, kappa, v_val,
			                         E_G, E0, L0, f, k_d, z_w, b_w, cw,
			                         dt = dt)
			utils::tail(traj$L, 1)
		}, numeric(1))

		# Normalise to control (C=0) for this draw
		L_ctrl <- L_final[1]  # c_seq starts at 0
		if (!is.finite(L_ctrl) || L_ctrl < 1e-12) L_ctrl <- NA_real_
		data.frame(
			concentration   = c_seq,
			relative        = L_final / L_ctrl,
			draw            = i
		)
	})
	pred_df <- do.call(rbind, pred_list)

	ggplot2::ggplot() +
		ggplot2::geom_line(
			data = pred_df,
			ggplot2::aes(x = .data$concentration, y = .data$relative,
			             group = .data$draw),
			alpha = 0.1, colour = "steelblue"
		) +
		ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
		                    colour = "grey50") +
		ggplot2::geom_point(
			data = obs_df,
			ggplot2::aes(x = .data$concentration, y = .data$relative),
			size = 3, colour = "red"
		) +
		ggplot2::theme_bw() +
		ggplot2::labs(
			title = "DEBtox Dose-Response (full model prediction)",
			x     = "Concentration",
			y     = expression(L[final] / L[control])
		) +
		ggplot2::scale_y_continuous(limits = c(0, 1.15))
}

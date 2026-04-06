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
#' \dontrun{
#' conc <- c("ctrl" = 0, "low" = 5, "mid" = 20, "high" = 100)
#' dat <- bdeb_data(growth = growth_df, concentration = conc)
#' mod <- bdeb_tox(dat, stress = "assimilation")
#' }
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
#' @return A named list with:
#'   - `draws`: posterior draws of EC50
#'   - `summary`: mean, median, sd, lower, upper
#'   - `NEC`: posterior summary of the no-effect concentration
#' @export
bdeb_ec50 <- function(fit, prob = 0.90) {
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
		NEC     = nec_summary
	)

	cli::cli_h3("DEBtox Effect Concentrations")
	print(result$summary, row.names = FALSE, digits = 3)

	invisible(result)
}

#' Plot DEBtox Dose-Response
#'
#' Creates a dose-response plot showing predicted growth inhibition
#' across concentrations with posterior uncertainty bands.
#'
#' @param fit A [bdeb_fit()] object from a DEBtox model.
#' @param endpoint Which endpoint to plot: `"growth"` or `"reproduction"`.
#'   Default `"growth"`.
#' @param n_draws Number of posterior draws. Default 100.
#' @return A ggplot2 object.
#' @export
plot_dose_response <- function(fit, endpoint = "growth", n_draws = 100) {
	if (!inherits(fit, "bdeb_fit") || fit$model$type != "debtox") {
		cli::cli_abort("Requires a fitted DEBtox model.")
	}

	draws <- posterior::as_draws_df(fit$fit$draws())
	n_total <- nrow(draws)
	idx <- sort(sample.int(n_total, min(n_draws, n_total)))

	# Get concentration groups
	C_w <- fit$model$stan_data$C_w
	N_groups <- length(C_w)

	# For each draw, compute predicted final length per concentration
	# This uses the last observation time's L_hat
	L_obs_mat <- fit$model$stan_data$L_obs
	N_obs <- fit$model$stan_data$N_obs

	# Collect observed endpoints (mean per group of last observation)
	obs_final <- vapply(seq_len(N_groups), function(g) {
		vals <- L_obs_mat[g, seq_len(N_obs[g])]
		vals <- vals[!is.nan(vals)]
		if (length(vals) > 0) utils::tail(vals, 1) else NA_real_
	}, numeric(1))

	obs_df <- data.frame(concentration = C_w, length = obs_final)

	# Predicted: use EC50 model for continuous curve
	c_seq <- seq(0, max(C_w) * 1.2, length.out = 100)

	pred_list <- lapply(idx, function(i) {
		z_w <- draws$z_w[i]
		b_w <- draws$b_w[i]
		# Stress factor at each concentration (steady state)
		s <- b_w * pmax(c_seq - z_w, 0)
		# Relative effect on assimilation
		effect <- pmax(1 - s, 0)
		data.frame(
			concentration = c_seq,
			relative_effect = effect,
			draw = i
		)
	})
	pred_df <- do.call(rbind, pred_list)

	ggplot2::ggplot() +
		ggplot2::geom_line(
			data = pred_df,
			ggplot2::aes(x = .data$concentration, y = .data$relative_effect,
			             group = .data$draw),
			alpha = 0.1, colour = "steelblue"
		) +
		ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
		ggplot2::geom_point(
			data = obs_df,
			ggplot2::aes(x = .data$concentration,
			             y = .data$length / max(.data$length, na.rm = TRUE)),
			size = 3, colour = "red"
		) +
		ggplot2::theme_bw() +
		ggplot2::labs(
			title = "DEBtox Dose-Response",
			x = "Concentration",
			y = "Relative Effect (1 = no effect)"
		) +
		ggplot2::scale_y_continuous(limits = c(0, 1.1))
}

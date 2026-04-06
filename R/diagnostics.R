#' MCMC Convergence Diagnostics
#'
#' Reports a comprehensive set of NUTS/HMC diagnostics following the
#' recommendations of Vehtari et al. (2021):
#' \describe{
#'   \item{Divergent transitions}{Indicate that the numerical leapfrog
#'     integrator encountered regions of high curvature.  Even a single
#'     divergence can bias the posterior.  Remedy: increase `adapt_delta`.}
#'   \item{Treedepth saturation}{The NUTS trajectory hit the maximum
#'     allowed tree depth, meaning it could not find a U-turn.  Remedy:
#'     increase `max_treedepth`.}
#'   \item{E-BFMI}{Energy Bayesian Fraction of Missing Information.
#'     Values below 0.3 indicate that the momentum resampling is
#'     inefficient (Betancourt, 2016).}
#'   \item{\eqn{\hat{R}}}{Split-\eqn{\hat{R}} convergence diagnostic.
#'     Values > 1.01 indicate incomplete mixing across chains.}
#'   \item{Bulk and tail ESS}{Effective sample size for the bulk and
#'     tails of the posterior.  Values below 400 suggest that posterior
#'     summaries may be unreliable.}
#' }
#'
#' @param fit A [bdeb_fit()] object.
#' @param pars Character vector of parameter names to report.  Default:
#'   all model parameters (excluding generated quantities such as
#'   `log_lik`, `L_rep`, and `lp__`).
#' @return Invisibly returns a list with components `n_divergent`,
#'   `n_max_treedepth`, `ebfmi`, and `summary` (a
#'   [posterior::summarise_draws()] tibble).
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Bürkner, P.-C. (2021). Rank-normalization, folding, and localization:
#' an improved \eqn{\hat{R}} for assessing convergence of MCMC.
#' *Bayesian Analysis*, 16(2), 667--718. \doi{10.1214/20-BA1221}
#'
#' Betancourt, M. (2016). Diagnosing biased inference with divergences.
#' Stan case study. \url{https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html}
#'
#' @export
bdeb_diagnose <- function(fit, pars = NULL) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	cli::cli_h2("BDEB Diagnostics")

	# --- CmdStan diagnostics ---
	diag <- fit$fit$diagnostic_summary(quiet = TRUE)

	n_div  <- sum(diag$num_divergent)
	n_tree <- sum(diag$num_max_treedepth)
	ebfmi  <- diag$ebfmi

	if (n_div > 0) {
		cli::cli_alert_danger("Divergent transitions: {n_div}")
		cli::cli_alert_info("Consider: increase adapt_delta, reparameterise, or tighten priors.")
	} else {
		cli::cli_alert_success("No divergent transitions.")
	}

	if (n_tree > 0) {
		cli::cli_alert_warning("Max treedepth saturated: {n_tree} times.")
		cli::cli_alert_info("Consider: increase max_treedepth.")
	} else {
		cli::cli_alert_success("Treedepth OK.")
	}

	low_ebfmi <- which(ebfmi < 0.3)
	if (length(low_ebfmi) > 0) {
		cli::cli_alert_warning("Low E-BFMI for chain(s): {low_ebfmi}")
	} else {
		cli::cli_alert_success("E-BFMI OK (all > 0.3).")
	}

	# --- Parameter-level diagnostics ---
	draws <- posterior::as_draws_df(fit$fit$draws())

	if (is.null(pars)) {
		# Get model parameter names (exclude log_lik, *_rep, lp__)
		all_vars <- posterior::variables(draws)
		pars <- all_vars[!grepl("^(log_lik|L_rep|R_rep|lp__|p_Am_new)", all_vars)]
	}

	summ <- posterior::summarise_draws(
		posterior::subset_draws(draws, variable = pars),
		"mean", "sd", "median",
		"q5" = ~ quantile(.x, 0.05),
		"q95" = ~ quantile(.x, 0.95),
		"rhat",
		"ess_bulk",
		"ess_tail"
	)

	cli::cli_h3("Parameter Summary")

	# Check for problematic Rhat
	bad_rhat <- summ$variable[!is.na(summ$rhat) & summ$rhat > 1.01]
	if (length(bad_rhat) > 0) {
		cli::cli_alert_danger("R-hat > 1.01 for: {paste(bad_rhat, collapse = ', ')}")
	} else {
		cli::cli_alert_success("All R-hat < 1.01.")
	}

	# Check for low ESS
	low_ess <- summ$variable[!is.na(summ$ess_bulk) & summ$ess_bulk < 400]
	if (length(low_ess) > 0) {
		cli::cli_alert_warning("Low bulk ESS (<400) for: {paste(low_ess, collapse = ', ')}")
	} else {
		cli::cli_alert_success("Bulk ESS adequate (>400) for all parameters.")
	}

	print(as.data.frame(summ), digits = 3, row.names = FALSE)

	invisible(list(
		n_divergent    = n_div,
		n_max_treedepth = n_tree,
		ebfmi          = ebfmi,
		summary        = summ
	))
}

#' Posterior Summary for BDEB Parameters
#'
#' Returns a tidy summary table of posterior draws for model parameters
#' and optionally derived quantities.
#'
#' @param fit A [bdeb_fit()] object.
#' @param pars Character vector of parameter names. Default: all model
#'   parameters.
#' @param prob Probability for credible intervals. Default 0.90 (5th/95th
#'   percentiles).
#' @param ... Ignored.
#' @return A `posterior::draws_summary` data frame.
#' @export
bdeb_summary <- function(fit, pars = NULL, prob = 0.90, ...) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	draws <- posterior::as_draws_df(fit$fit$draws())

	if (is.null(pars)) {
		all_vars <- posterior::variables(draws)
		pars <- all_vars[!grepl("^(log_lik|lp__)", all_vars)]
	}

	alpha <- (1 - prob) / 2
	posterior::summarise_draws(
		posterior::subset_draws(draws, variable = pars),
		"mean", "sd", "median",
		"lower" = ~ quantile(.x, alpha),
		"upper" = ~ quantile(.x, 1 - alpha),
		"rhat",
		"ess_bulk",
		"ess_tail"
	)
}

#' Compute Derived Biological Quantities from the Posterior
#'
#' Transforms the raw DEB parameter draws into biologically interpretable
#' quantities, automatically propagating parameter uncertainty.  The
#' formulas follow Kooijman (2010, Ch. 3) and the AmP parameter
#' definitions (Marques et al., 2018):
#'
#' \describe{
#'   \item{`"L_inf"`}{Ultimate structural length at food level \eqn{f}:
#'     \deqn{L_\infty = f \, \kappa \, \{p_{Am}\} / [p_M]}
#'     This is the asymptotic body length when \eqn{dV/dt = 0}.}
#'   \item{`"k_M"`}{Somatic maintenance rate constant:
#'     \deqn{k_M = [p_M] / [E_G]}
#'     The ratio of maintenance costs to structural costs (units: d\eqn{^{-1}}).}
#'   \item{`"growth_rate"`}{Von Bertalanffy growth rate:
#'     \deqn{\dot{r}_B = \frac{v}{3} \cdot \frac{[p_M]}{\kappa \, [E_G]}}
#'     Describes how quickly the organism approaches \eqn{L_\infty}.}
#' }
#'
#' @param fit A [bdeb_fit()] object.
#' @param quantities Character vector of quantities to compute.  One or
#'   more of `"L_inf"`, `"k_M"`, `"growth_rate"`.
#' @param f Scaled functional response \eqn{f \in (0,1]} for computing
#'   food-dependent quantities.  Default 1 (ad libitum).
#' @return A [posterior::draws_df] with one column per requested quantity
#'   and one row per posterior draw.
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press.
#' \doi{10.1017/CBO9780511805400}
#'
#' Marques, G.M., Augustine, S., Lika, K., Pecquerie, L., Domingos, T.
#' and Kooijman, S.A.L.M. (2018). The AmP project: comparing species on
#' the basis of dynamic energy budget parameters. *PLOS Computational
#' Biology*, 14(5), e1006100. \doi{10.1371/journal.pcbi.1006100}
#' @export
bdeb_derived <- function(fit,
                         quantities = c("L_inf", "k_M", "growth_rate"),
                         f = 1.0) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	draws <- posterior::as_draws_df(fit$fit$draws())

	result <- data.frame(.draw = seq_len(nrow(draws)))

	# Extract parameter draws
	get_par <- function(name) {
		if (name %in% names(draws)) return(draws[[name]])
		NULL
	}

	p_Am  <- get_par("p_Am")
	p_M   <- get_par("p_M")
	kappa <- get_par("kappa")
	v     <- get_par("v")
	E_G   <- get_par("E_G")

	if (is.null(p_Am)) {
		# Hierarchical: use population mean
		mu <- get_par("mu_log_p_Am")
		if (!is.null(mu)) p_Am <- exp(mu)
	}

	if ("L_inf" %in% quantities && !is.null(p_Am) && !is.null(p_M)) {
		# L_inf = f * kappa * p_Am / p_M  (from von Bertalanffy at steady state)
		# More precisely: L_inf = (f * p_Am * kappa) / p_M
		# in volumetric terms then L_inf = (f * kappa * {p_Am} / [p_M])^(1/3)
		# Simplified for direct parameters:
		result$L_inf <- (f * kappa * p_Am / p_M)
	}

	if ("k_M" %in% quantities && !is.null(p_M) && !is.null(E_G)) {
		result$k_M <- p_M / E_G
	}

	if ("growth_rate" %in% quantities && !is.null(v) && !is.null(p_M) &&
	    !is.null(kappa) && !is.null(E_G)) {
		# Von Bertalanffy rate: k_vB = v / 3 * p_M / (kappa * E_G)
		result$growth_rate <- (v / 3) * (p_M / (kappa * E_G))
	}

	posterior::as_draws_df(result)
}

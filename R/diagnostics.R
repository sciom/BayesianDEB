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
#' @param verbose Logical; if `TRUE` (default) the diagnostic messages
#'   and parameter summary table are printed via `cli` functions and
#'   [message()].  All output is suppressible with
#'   [suppressMessages()].  Set to `FALSE` for a silent run (the
#'   invisible return value is unchanged).
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
bdeb_diagnose <- function(fit, pars = NULL, verbose = TRUE) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	if (verbose) cli::cli_h2("BDEB Diagnostics")

	# --- CmdStan diagnostics ---
	diag <- fit$fit$diagnostic_summary(quiet = TRUE)

	n_div  <- sum(diag$num_divergent)
	n_tree <- sum(diag$num_max_treedepth)
	ebfmi  <- diag$ebfmi

	if (verbose) {
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

	if (verbose) {
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

		# Route the summary table through message() so it can be
		# silenced with suppressMessages() — CRAN requirement.
		tbl_lines <- utils::capture.output(
			print(as.data.frame(summ), digits = 3, row.names = FALSE)
		)
		cli::cli_verbatim(tbl_lines)
	}

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
		pars <- all_vars[!grepl("^(log_lik|L_hat|L_rep|R_hat|R_rep|lp__|p_Am_new)", all_vars)]
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

#' LOO Cross-Validation for Model Comparison
#'
#' Computes approximate leave-one-out cross-validation (LOO-CV) using
#' Pareto-smoothed importance sampling (PSIS; Vehtari et al., 2017).
#' Requires that the Stan model includes a `log_lik` array in the
#' `generated quantities` block.
#'
#' Currently supported for `"individual"` and `"growth_repro"` models
#' only.  The `"hierarchical"` and `"debtox"` models do not store
#' per-observation `log_lik` in `generated quantities`:
#' `"hierarchical"` computes the likelihood inside `reduce_sum`
#' (no access to individual contributions outside the function);
#' `"debtox"` uses the same `reduce_sum` approach and its
#' `generated quantities` block only computes EC50/NEC.
#' Adding per-group `log_lik` to DEBtox is planned for a future
#' version.  An informative error is raised for unsupported types.
#'
#' @section Conditional independence assumption:
#' For `"growth_repro"` models with `endpoint = "all"`, growth and
#' reproduction observations are concatenated and each treated as an
#' independent data point.  This is valid because the two endpoints are
#' **conditionally independent given the latent DEB process** (growth
#' observations depend only on \eqn{V(t)}, reproduction counts depend
#' only on \eqn{\Delta E_R}; given the ODE solution, they share no
#' additional error).  Use `endpoint = "growth"` or
#' `endpoint = "reproduction"` to compute LOO for a single endpoint.
#'
#' @param fit A [bdeb_fit()] object.
#' @param endpoint Which log-likelihood to use for `"growth_repro"`
#'   models: `"all"` (default, concatenates growth + reproduction),
#'   `"growth"`, or `"reproduction"`.  Ignored for `"individual"` models.
#' @param ... Additional arguments passed to [loo::loo()].
#' @return A `loo` object (see [loo::loo()]).
#'
#' @references
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian
#' model evaluation using leave-one-out cross-validation and WAIC.
#' *Statistics and Computing*, 27(5), 1413--1432.
#' \doi{10.1007/s11222-016-9696-4}
#' @export
bdeb_loo <- function(fit, endpoint = c("all", "growth", "reproduction"), ...) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	if (!requireNamespace("loo", quietly = TRUE)) {
		cli::cli_abort(c(
			"x" = "{.pkg loo} is required for LOO-CV.",
			"i" = 'Install with: install.packages("loo")'
		))
	}

	mtype <- fit$model$type

	if (mtype %in% c("hierarchical", "debtox")) {
		cli::cli_abort(c(
			"x" = "LOO-CV is not available for {.val {mtype}} models.",
			"i" = "These models use {.code reduce_sum} and do not store per-observation {.var log_lik}.",
			"i" = "Use {.fun bdeb_diagnose} and posterior predictive checks for model assessment."
		))
	}

	endpoint <- match.arg(endpoint)
	draws <- fit$fit$draws(format = "draws_matrix")

	# Find log_lik columns based on endpoint selection
	if (mtype == "individual") {
		ll_vars <- grep("^log_lik\\[", colnames(draws), value = TRUE)
	} else {
		# growth_repro has log_lik_L and log_lik_R
		ll_L <- grep("^log_lik_L", colnames(draws), value = TRUE)
		ll_R <- grep("^log_lik_R", colnames(draws), value = TRUE)

		ll_vars <- switch(endpoint,
			all          = c(ll_L, ll_R),
			growth       = ll_L,
			reproduction = ll_R
		)
	}

	if (length(ll_vars) == 0) {
		cli::cli_abort("No {.var log_lik} variables found for endpoint {.val {endpoint}}.")
	}

	log_lik <- draws[, ll_vars, drop = FALSE]
	loo::loo(log_lik, ...)
}

#' Compute Derived Biological Quantities from the Posterior
#'
#' Transforms the raw DEB parameter draws into biologically interpretable
#' quantities, automatically propagating parameter uncertainty.  All
#' formulas follow Kooijman (2010, Ch. 3, Table 3.1) strictly.
#'
#' @section DEB length terminology:
#' DEB theory distinguishes several length measures.  This function
#' returns **structural length** \eqn{L = V^{1/3}}, which is the cube
#' root of structural volume.  Structural length is *not* the same as
#' **physical (observed) length** \eqn{L_w}, which relates to \eqn{L}
#' via the shape coefficient: \eqn{L_w = L / \delta_M}.  The shape
#' coefficient \eqn{\delta_M} is species-specific (typically 0.1--0.5)
#' and is not estimated by this package.  If your data are physical
#' lengths, you should either (a) convert to structural length before
#' fitting, or (b) divide `L_inf` by your species' \eqn{\delta_M} to
#' obtain the physical asymptotic length.
#'
#' @section Available quantities:
#' \describe{
#'   \item{`"L_m"`}{Maximum structural length at \eqn{f = 1}
#'     (Kooijman, 2010, Table 3.1):
#'     \deqn{L_m = \kappa \, \{p_{Am}\} / [p_M]}
#'     This is a **species-level constant** independent of food.
#'     Units: cm.}
#'   \item{`"L_inf"`}{Ultimate structural length at food level \eqn{f}
#'     (Kooijman, 2010, Eq. 3.4):
#'     \deqn{L_i = f \cdot L_m = f \, \kappa \, \{p_{Am}\} / [p_M]}
#'     The asymptotic length when \eqn{dV/dt = 0} at constant food.
#'     Depends on \eqn{f}.
#'     Units: cm.  Dimensional check:
#'     \eqn{(-)(-)(\text{J d}^{-1}\text{cm}^{-2}) /
#'     (\text{J d}^{-1}\text{cm}^{-3}) = \text{cm}}.}
#'   \item{`"k_M"`}{Somatic maintenance rate constant:
#'     \deqn{k_M = [p_M] / [E_G]}
#'     Units: d\eqn{^{-1}}.}
#'   \item{`"growth_rate"`}{Von Bertalanffy growth rate
#'     (Kooijman, 2010, Eq. 3.23):
#'     \deqn{\dot{r}_B = \frac{k_M \, g}{3\,(f + g)}}
#'     where \eqn{g = [E_G] \, v / (\kappa \, \{p_{Am}\})} is the energy
#'     investment ratio.  Depends on \eqn{f}.  Units: d\eqn{^{-1}}.}
#'   \item{`"g"`}{Energy investment ratio (Kooijman, 2010, Table 3.1):
#'     \deqn{g = [E_G] \, v / (\kappa \, \{p_{Am}\})}
#'     Dimensionless.  Large \eqn{g} means growth is expensive relative
#'     to reserve turnover.}
#' }
#'
#' @section Reference example:
#' For *Eisenia fetida* with AmP parameters \eqn{\{p_{Am}\} = 5.0},
#' \eqn{[p_M] = 0.5}, \eqn{\kappa = 0.75}:
#' \eqn{L_m = 0.75 \times 5.0 / 0.5 = 7.5} cm (structural).
#' With \eqn{\delta_M \approx 0.24}, the physical maximum length would
#' be \eqn{L_m / \delta_M \approx 31} mm, consistent with observations.
#'
#' @param fit A [bdeb_fit()] object.
#' @param quantities Character vector of quantities to compute.  One or
#'   more of `"L_m"`, `"L_inf"`, `"k_M"`, `"growth_rate"`, `"g"`.
#' @param f Scaled functional response \eqn{f \in (0,1]} for computing
#'   food-dependent quantities (`"L_inf"`, `"growth_rate"`).
#'   Default 1 (ad libitum).
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

	if ("L_m" %in% quantities && !is.null(p_Am) && !is.null(p_M)) {
		# Maximum structural length (at f=1): L_m = kappa * {p_Am} / [p_M]
		# Kooijman (2010), Table 3.1.  Units: cm.
		result$L_m <- kappa * p_Am / p_M
	}

	if ("L_inf" %in% quantities && !is.null(p_Am) && !is.null(p_M)) {
		# Ultimate structural length at food level f: L_i = f * L_m
		# Kooijman (2010), Eq. 3.4.  Units: cm.
		# NOTE: this is structural length (V^{1/3}), not physical length.
		# Physical length = L / delta_M (shape coefficient, species-specific).
		result$L_inf <- f * kappa * p_Am / p_M
	}

	if ("k_M" %in% quantities && !is.null(p_M) && !is.null(E_G)) {
		# Somatic maintenance rate constant: k_M = [p_M] / [E_G]
		# Units: d^-1
		result$k_M <- p_M / E_G
	}

	if ("g" %in% quantities && !is.null(v) && !is.null(E_G) &&
	    !is.null(kappa) && !is.null(p_Am)) {
		# Energy investment ratio: g = [E_G] * v / (kappa * {p_Am})
		# Kooijman (2010), Table 3.1.  Dimensionless.
		result$g <- E_G * v / (kappa * p_Am)
	}

	if ("growth_rate" %in% quantities && !is.null(v) && !is.null(p_M) &&
	    !is.null(kappa) && !is.null(E_G) && !is.null(p_Am)) {
		# Von Bertalanffy growth rate: r_B = k_M * g / (3 * (f + g))
		# Kooijman (2010), Eq. 3.23.  Units: d^-1.
		# where g = [E_G]*v / (kappa*{p_Am}), k_M = [p_M]/[E_G]
		g <- E_G * v / (kappa * p_Am)
		k_M <- p_M / E_G
		result$growth_rate <- k_M * g / (3 * (f + g))
	}

	posterior::as_draws_df(result)
}

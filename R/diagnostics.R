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
#' Returns a `bdeb_diagnostics` S3 object.  The object has dedicated
#' [print()][print.bdeb_diagnostics],
#' [summary()][summary.bdeb_diagnostics] and
#' [plot()][plot.bdeb_diagnostics] methods.  When called interactively
#' the print method is invoked automatically; assign the result
#' (`d <- bdeb_diagnose(fit)`) to suppress output.
#'
#' @param fit A [bdeb_fit()] object.
#' @param pars Character vector of parameter names to report.  Default:
#'   all model parameters (excluding generated quantities such as
#'   `log_lik`, `L_rep`, and `lp__`).
#' @return An object of class `bdeb_diagnostics` with components
#'   `n_divergent`, `n_max_treedepth`, `ebfmi`, `summary` (a
#'   [posterior::summarise_draws()] tibble), `pars`, and `model_type`.
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
#' @seealso [print.bdeb_diagnostics()], [summary.bdeb_diagnostics()],
#'   [plot.bdeb_diagnostics()]
#' @export
#' @examples
#' # Requires the CmdStan toolchain (Suggests: cmdstanr); gated on its
#' # availability and wrapped in \donttest{} so example("bdeb_diagnose")
#' # runs it when a toolchain is present, mirroring bdeb_fit().
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(),
#'                     error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   d <- bdeb_diagnose(fit)
#'   print(d)
#'   summary(d)
#'   plot(d, type = "rhat")
#' }
#' }
bdeb_diagnose <- function(fit, pars = NULL) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}
	if (!is.null(pars) && !is.character(pars)) {
		cli::cli_abort("{.arg pars} must be a character vector or NULL.")
	}

	algo <- if (is.null(fit$algorithm)) "sampling" else fit$algorithm
	if (!identical(algo, "sampling")) {
		cli::cli_abort(c(
			"x" = "{.fn bdeb_diagnose} requires a sampling (NUTS) fit.",
			"i" = "This fit was produced with {.code algorithm = {.val {algo}}}, which does not provide R-hat, divergent transitions, or treedepth diagnostics.",
			"i" = "Re-fit with {.code algorithm = \"sampling\"} for full diagnostics."
		))
	}

	diag <- fit$fit$diagnostic_summary(quiet = TRUE)

	draws <- posterior::as_draws_df(fit$fit$draws())

	if (is.null(pars)) {
		all_vars <- posterior::variables(draws)
		pars <- all_vars[!grepl("^(log_lik|L_rep|R_rep|lp__|p_Am_new)", all_vars)]
	}

	summ <- posterior::summarise_draws(
		posterior::subset_draws(draws, variable = pars),
		"mean", "sd", "median",
		"q5" = ~ quantile(.x, 0.05, na.rm = TRUE),
		"q95" = ~ quantile(.x, 0.95, na.rm = TRUE),
		"rhat",
		"ess_bulk",
		"ess_tail"
	)

	out <- list(
		n_divergent     = sum(diag$num_divergent),
		n_max_treedepth = sum(diag$num_max_treedepth),
		ebfmi           = diag$ebfmi,
		summary         = summ,
		pars            = pars,
		model_type      = fit$model$type
	)
	structure(out, class = "bdeb_diagnostics")
}

#' Print a BDEB Diagnostics Report
#'
#' Default printing for [bdeb_diagnose()] output.  Displays divergence /
#' treedepth / E-BFMI alerts, R-hat and ESS warnings, and a compact
#' parameter summary table.  To keep the on-screen output short, the
#' per-time-point latent states (`x_sol[i,j]`, `L_hat[i]`, ...) are hidden
#' by default; the scalar model parameters are always shown.  Output uses
#' [cli] alerts and is therefore silenceable via [cli::cli_inform()] sinks.
#'
#' @param x A `bdeb_diagnostics` object.
#' @param full Logical; if `TRUE`, also print the latent-state rows that
#'   are hidden by default.  The complete table is always available via
#'   `summary(x)$table`.
#' @param ... Unused.
#' @return The input object, invisibly.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   print(bdeb_diagnose(fit))
#' }
#' }
print.bdeb_diagnostics <- function(x, full = FALSE, ...) {
	cli::cli_h2("BDEB Diagnostics ({x$model_type})")

	if (x$n_divergent > 0) {
		cli::cli_alert_danger("Divergent transitions: {x$n_divergent}")
		cli::cli_alert_info("Consider: increase adapt_delta, reparameterise, or tighten priors.")
	} else {
		cli::cli_alert_success("No divergent transitions.")
	}

	if (x$n_max_treedepth > 0) {
		cli::cli_alert_warning("Max treedepth saturated: {x$n_max_treedepth} times.")
		cli::cli_alert_info("Consider: increase max_treedepth.")
	} else {
		cli::cli_alert_success("Treedepth OK.")
	}

	low_ebfmi <- which(x$ebfmi < 0.3)
	if (length(low_ebfmi) > 0) {
		cli::cli_alert_warning("Low E-BFMI for chain(s): {low_ebfmi}")
	} else {
		cli::cli_alert_success("E-BFMI OK (all > 0.3).")
	}

	cli::cli_h3("Parameter Summary")

	bad_rhat <- x$summary$variable[!is.na(x$summary$rhat) & x$summary$rhat > 1.01]
	if (length(bad_rhat) > 0) {
		cli::cli_alert_danger("R-hat > 1.01 for: {paste(bad_rhat, collapse = ', ')}")
	} else {
		cli::cli_alert_success("All R-hat < 1.01.")
	}

	low_ess <- x$summary$variable[!is.na(x$summary$ess_bulk) & x$summary$ess_bulk < 400]
	if (length(low_ess) > 0) {
		cli::cli_alert_warning("Low bulk ESS (<400) for: {paste(low_ess, collapse = ', ')}")
	} else {
		cli::cli_alert_success("Bulk ESS adequate (>400) for all parameters.")
	}

	# Keep the on-screen table short: by default show only the scalar
	# model parameters and hide the per-time-point latent states
	# (x_sol[i,j], L_hat[i], ...), which can run to dozens of rows.
	tbl <- as.data.frame(x$summary)
	is_latent <- grepl("\\[", tbl$variable)
	n_hidden  <- sum(is_latent)
	show_tbl  <- if (full || n_hidden == 0L) tbl else tbl[!is_latent, , drop = FALSE]

	tbl_lines <- utils::capture.output(
		print(show_tbl, digits = 3, row.names = FALSE)
	)
	cli::cli_verbatim(tbl_lines)

	if (!full && n_hidden > 0L) {
		cli::cli_alert_info(
			"{n_hidden} latent-state row{?s} hidden; use {.code print(x, full = TRUE)} \\
			 or {.code summary(x)$table} to see all.")
	}

	invisible(x)
}

#' Compact Summary of a BDEB Diagnostics Report
#'
#' Returns counts of problematic parameters (divergences, treedepth
#' saturations, low-EBFMI chains, R-hat > 1.01, ESS-bulk < 400) suitable
#' for a one-line health check or programmatic gating.
#'
#' @param object A `bdeb_diagnostics` object.
#' @param ... Unused.
#' @return An object of class `summary.bdeb_diagnostics` (a list).
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   summary(bdeb_diagnose(fit))
#' }
#' }
summary.bdeb_diagnostics <- function(object, ...) {
	bad_rhat <- object$summary$variable[
		!is.na(object$summary$rhat) & object$summary$rhat > 1.01]
	low_ess <- object$summary$variable[
		!is.na(object$summary$ess_bulk) & object$summary$ess_bulk < 400]
	low_ebfmi <- which(object$ebfmi < 0.3)

	out <- list(
		model_type      = object$model_type,
		n_pars          = length(object$pars),
		n_divergent     = object$n_divergent,
		n_max_treedepth = object$n_max_treedepth,
		n_low_ebfmi     = length(low_ebfmi),
		n_bad_rhat      = length(bad_rhat),
		n_low_ess       = length(low_ess),
		bad_rhat        = bad_rhat,
		low_ess         = low_ess,
		table           = object$summary
	)
	structure(out, class = "summary.bdeb_diagnostics")
}

#' @export
print.summary.bdeb_diagnostics <- function(x, ...) {
	cli::cli_h2("BDEB Diagnostics summary ({x$model_type})")
	cli::cli_li("Parameters monitored: {.val {x$n_pars}}")
	cli::cli_li("Divergent transitions: {.val {x$n_divergent}}")
	cli::cli_li("Treedepth saturations: {.val {x$n_max_treedepth}}")
	cli::cli_li("Chains with low E-BFMI: {.val {x$n_low_ebfmi}}")
	cli::cli_li("Parameters with R-hat > 1.01: {.val {x$n_bad_rhat}}")
	cli::cli_li("Parameters with ESS-bulk < 400: {.val {x$n_low_ess}}")
	invisible(x)
}

#' Plot Convergence Diagnostics
#'
#' Visualises the per-parameter R-hat or ESS-bulk values from a
#' [bdeb_diagnose()] object.  A dashed red reference line is drawn at the
#' Vehtari et al. (2021) threshold (R-hat = 1.01, ESS-bulk = 400).
#'
#' To keep the plot short and readable, the per-time-point latent states
#' (`x_sol[i,j]`, `L_hat[i]`, ...) are hidden by default and only the
#' scalar model parameters are shown, mirroring [print.bdeb_diagnostics()];
#' set `full = TRUE` to plot every monitored quantity.
#'
#' @param x A `bdeb_diagnostics` object.
#' @param type One of `"rhat"` (default) or `"ess"`.
#' @param full Logical; if `TRUE`, also plot the latent-state rows that
#'   are hidden by default.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   plot(bdeb_diagnose(fit), type = "ess")
#' }
#' }
plot.bdeb_diagnostics <- function(x, type = c("rhat", "ess"), full = FALSE, ...) {
	type <- match.arg(type)
	df <- as.data.frame(x$summary)

	# Keep the plot short: by default drop the per-time-point latent
	# states (x_sol[i,j], L_hat[i], ...) and show only scalar model
	# parameters, mirroring print.bdeb_diagnostics().
	is_latent <- grepl("\\[", df$variable)
	n_hidden  <- sum(is_latent)
	if (!full && n_hidden > 0L) {
		df <- df[!is_latent, , drop = FALSE]
	}
	hidden_note <- if (!full && n_hidden > 0L) {
		sprintf("  (%d latent-state row%s hidden; full = TRUE to show all)",
		        n_hidden, if (n_hidden == 1L) "" else "s")
	} else ""

	if (type == "rhat") {
		ggplot2::ggplot(
			df,
			ggplot2::aes(x = .data$rhat,
			             y = stats::reorder(.data$variable, .data$rhat))
		) +
			ggplot2::geom_point(size = 2) +
			ggplot2::geom_vline(xintercept = 1.01,
			                     linetype = "dashed", colour = "red") +
			ggplot2::labs(
				x = expression(hat(R)), y = NULL,
				title = "Convergence: split-Rhat",
				subtitle = paste0("Dashed line at 1.01 (Vehtari et al., 2021)",
				                  hidden_note)
			) +
			ggplot2::theme_minimal()
	} else {
		ggplot2::ggplot(
			df,
			ggplot2::aes(x = .data$ess_bulk,
			             y = stats::reorder(.data$variable, .data$ess_bulk))
		) +
			ggplot2::geom_point(size = 2) +
			ggplot2::geom_vline(xintercept = 400,
			                     linetype = "dashed", colour = "red") +
			ggplot2::labs(
				x = "ESS-bulk", y = NULL,
				title = "Effective sample size (bulk)",
				subtitle = paste0("Dashed line at 400", hidden_note)
			) +
			ggplot2::theme_minimal()
	}
}

#' Posterior Summary for BDEB Parameters (deprecated)
#'
#' @description
#' `bdeb_summary()` is deprecated as of BayesianDEB 0.2.0.  Call
#' [summary()] on a [bdeb_fit()] object instead; the two are
#' equivalent.  This wrapper will be removed in a future release.
#'
#' @param fit A [bdeb_fit()] object.
#' @param pars Character vector of parameter names.  Forwarded to
#'   the [summary()] method on `bdeb_fit`.
#' @param prob Probability for the central credible interval.  Forwarded
#'   to the [summary()] method on `bdeb_fit`.
#' @param ... Ignored.
#' @return A `posterior::draws_summary` data frame.
#' @keywords internal
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   bdeb_summary(fit)  # equivalent to summary(fit); deprecated
#' }
#' }
bdeb_summary <- function(fit, pars = NULL, prob = 0.90, ...) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}
	.Deprecated(
		new = "summary",
		package = "BayesianDEB",
		msg = "`bdeb_summary()` is deprecated; use `summary(fit)` instead."
	)
	summary(fit, pars = pars, prob = prob, ...)
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
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   bdeb_loo(fit)
#' }
#' }
bdeb_loo <- function(fit, endpoint = c("all", "growth", "reproduction"), ...) {
	if (!inherits(fit, "bdeb_fit")) {
		cli::cli_abort("{.arg fit} must be a {.cls bdeb_fit} object.")
	}

	algo <- if (is.null(fit$algorithm)) "sampling" else fit$algorithm
	if (!identical(algo, "sampling")) {
		cli::cli_abort(c(
			"x" = "{.fn bdeb_loo} requires a sampling (NUTS) fit; PSIS-LOO is not defined for variational approximations.",
			"i" = "Re-fit with {.code algorithm = \"sampling\"} for cross-validation."
		))
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
#' @param object A [bdeb_fit()] object.
#' @param ... Additional arguments passed to methods.
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
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
#'   data(eisenia_growth)
#'   dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#'   fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
#'                   chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                   refresh = 0)
#'   bdeb_derived(fit, quantities = c("L_inf", "k_M"))
#' }
#' }
bdeb_derived <- function(object, ...) {
	UseMethod("bdeb_derived")
}

#' @rdname bdeb_derived
#' @export
bdeb_derived.default <- function(object, ...) {
	cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
}

#' @rdname bdeb_derived
#' @export
bdeb_derived.bdeb_fit <- function(object,
                                  quantities = c("L_inf", "k_M", "growth_rate"),
                                  f = 1.0,
                                  ...) {
	draws <- posterior::as_draws_df(object$fit$draws())

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

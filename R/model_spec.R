#' Specify a BDEB Model
#'
#' Creates a model specification that binds together the prepared data, the
#' DEB process model (encoded as a pre-written Stan program), the prior
#' distributions, and the observation model.  Internally it calls the
#' appropriate `build_stan_data_*()` function to assemble the list of
#' named values that Stan expects.  The Stan program is *not* compiled at
#' this stage — compilation and sampling happen in [bdeb_fit()].
#'
#' The four model types correspond to increasingly complex DEB formulations:
#'
#' \describe{
#'   \item{`"individual"`}{Standard 2-state DEB (reserve \eqn{E}, structure
#'     \eqn{V}), Kooijman (2010, Ch. 2).  The ODE is solved with Stan's
#'     `ode_bdf` (stiff BDF solver) at tolerances \eqn{10^{-6}}.}
#'   \item{`"growth_repro"`}{3-state model adding the reproduction buffer
#'     \eqn{E_R}.  Offspring counts are modelled as
#'     \eqn{R_i \sim \mathrm{NegBin}(k_R \Delta E_R, \phi)},
#'     where \eqn{\Delta E_R = E_R(t_{\mathrm{end}}) - E_R(t_{\mathrm{start}})}.}
#'   \item{`"hierarchical"`}{2-state model with a lognormal random effect on
#'     \eqn{\{p_{Am}\}}: \eqn{\log p_{Am,j} = \mu + \sigma z_j},
#'     \eqn{z_j \sim N(0,1)} (non-centred).  Shared parameters
#'     \eqn{[p_M], \kappa, v, [E_G]} are estimated from all individuals
#'     jointly (partial pooling; Gelman & Hill, 2006).}
#'   \item{`"debtox"`}{4-state TKTD model following Jager et al. (2006).
#'     Adds scaled damage \eqn{D_w} with
#'     \eqn{dD_w/dt = k_d(\max(C_w - z_w, 0) - D_w)}.
#'     The \eqn{\mathrm{EC}_{50}} is computed analytically as
#'     \eqn{z_w + 0.5/b_w}.}
#' }
#'
#' @param data A [bdeb_data()] object.
#' @param type Model type: `"individual"`, `"growth_repro"`,
#'   `"hierarchical"`, or `"debtox"`.
#' @param priors Named list of `bdeb_prior` objects (see [prior_lognormal()],
#'   [prior_beta()], etc.). Entries not supplied are filled from
#'   [prior_default()].
#' @param observation Named list of observation model specs for each endpoint.
#'   Default: `list(growth = obs_normal(), reproduction = obs_negbinom())`.
#' @param temperature Optional list with components `T` (observation
#'   temperature in Kelvin), `T_ref` (reference temperature), and
#'   `T_A` (Arrhenius temperature).  If provided, rate parameters
#'   are corrected by the factor \eqn{c_T = \exp(T_A/T_{\mathrm{ref}} - T_A/T)};
#'   see [arrhenius()].
#' @return A `bdeb_model` object (S3 list).
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press.
#' \doi{10.1017/CBO9780511805400}
#'
#' Gelman, A. and Hill, J. (2006). *Data Analysis Using Regression and
#' Multilevel/Hierarchical Models*. Cambridge University Press.
#'
#' Jager, T., Heugens, E.H.W. and Kooijman, S.A.L.M. (2006). Making
#' sense of ecotoxicological test results: towards application of
#' process-based models. *Ecotoxicology*, 15(3), 305--314.
#' \doi{10.1007/s10646-006-0060-x}
#'
#' @export
bdeb_model <- function(data,
                       type = c("individual", "growth_repro",
                                "hierarchical", "debtox"),
                       priors = list(),
                       observation = list(),
                       temperature = NULL) {

	if (!inherits(data, "bdeb_data")) {
		cli::cli_abort("{.arg data} must be a {.cls bdeb_data} object from {.fun bdeb_data}.")
	}

	type <- match.arg(type)

	# Validate type vs data
	if (type == "individual" && data$n_ind > 1) {
		cli::cli_alert_warning(
			"Individual model with {data$n_ind} individuals. Only first will be used."
		)
	}
	if (type == "growth_repro" && !("reproduction" %in% data$endpoints)) {
		cli::cli_abort("growth_repro model requires reproduction data in {.fun bdeb_data}.")
	}
	if (type == "debtox" && is.null(data$concentration)) {
		cli::cli_abort("debtox model requires {.arg concentration} in {.fun bdeb_data}.")
	}

	# Fill defaults
	defaults <- prior_default(type)
	for (nm in names(defaults)) {
		if (is.null(priors[[nm]])) priors[[nm]] <- defaults[[nm]]
	}

	# Observation models — defaults
	if (is.null(observation$growth))       observation$growth <- obs_normal()
	if (is.null(observation$reproduction)) observation$reproduction <- obs_negbinom()

	# Stan model selection
	stan_model_name <- switch(type,
		individual   = "bdeb_individual_growth",
		growth_repro = "bdeb_growth_repro",
		hierarchical = "bdeb_hierarchical_growth",
		debtox       = "bdeb_debtox"
	)

	# Build Stan data
	stan_data <- switch(type,
		individual   = build_stan_data_individual(data, priors),
		growth_repro = build_stan_data_growth_repro(data, priors),
		hierarchical = build_stan_data_hierarchical(data, priors),
		debtox       = build_stan_data_debtox(data, priors)
	)

	out <- list(
		data            = data,
		type            = type,
		priors          = priors,
		observation     = observation,
		temperature     = temperature,
		stan_model_name = stan_model_name,
		stan_data       = stan_data
	)

	structure(out, class = "bdeb_model")
}

#' @return The input object, invisibly.
#' @export
print.bdeb_model <- function(x, ...) {
	cli::cli_h2("BDEB Model Specification")
	cli::cli_alert_info("Type: {x$type}")
	cli::cli_alert_info("Stan model: {x$stan_model_name}")
	cli::cli_alert_info("Individuals: {x$data$n_ind}")
	cli::cli_alert_info("Endpoints: {paste(x$data$endpoints, collapse = ', ')}")

	cli::cli_h3("Priors")
	for (nm in names(x$priors)) {
		p <- x$priors[[nm]]
		desc <- switch(p$family,
			lognormal   = sprintf("LogNormal(%.1f, %.1f)", p$mu, p$sigma),
			normal      = sprintf("Normal(%.1f, %.1f)", p$mu, p$sigma),
			beta        = sprintf("Beta(%.1f, %.1f)", p$a, p$b),
			halfnormal  = sprintf("HalfNormal(%.2f)", p$sigma),
			halfcauchy  = sprintf("HalfCauchy(%.2f)", p$sigma),
			exponential = sprintf("Exponential(%.1f)", p$rate),
			"unknown"
		)
		cli::cli_alert("  {nm}: {desc}")
	}

	if (!is.null(x$temperature)) {
		cli::cli_h3("Temperature correction")
		cli::cli_alert("  T = {x$temperature$T} K, T_ref = {x$temperature$T_ref} K, T_A = {x$temperature$T_A} K")
	}

	invisible(x)
}

# --- Observation model specs ---

#' Observation Model Specifications
#'
#' @return A `bdeb_obs` object (list with class `"bdeb_obs"`).
#' @name observation_models
#' @examples
#' obs_normal()
#' obs_lognormal()
#' obs_negbinom()
NULL

#' @describeIn observation_models Gaussian observation error (default for growth)
#' @export
obs_normal <- function() {
	structure(list(family = "normal"), class = "bdeb_obs")
}

#' @describeIn observation_models Log-normal observation error (multiplicative)
#' @export
obs_lognormal <- function() {
	structure(list(family = "lognormal"), class = "bdeb_obs")
}

#' @describeIn observation_models Student-t observation error (robust to outliers)
#' @param nu Degrees of freedom. Default 5.
#' @export
obs_student_t <- function(nu = 5) {
	structure(list(family = "student_t", nu = nu), class = "bdeb_obs")
}

#' @describeIn observation_models Poisson observations (for count data)
#' @export
obs_poisson <- function() {
	structure(list(family = "poisson"), class = "bdeb_obs")
}

#' @describeIn observation_models Negative binomial observations (overdispersed counts)
#' @export
obs_negbinom <- function() {
	structure(list(family = "negbinom"), class = "bdeb_obs")
}

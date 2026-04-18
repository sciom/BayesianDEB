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
#'     jointly (partial pooling; Gelman & Hill, 2006).  Initial
#'     structural length \eqn{L_{0,j}} varies per individual; initial
#'     reserve density \eqn{E_0} is shared (assumes organisms start with
#'     the same reserve state, which is typical for lab-reared cohorts).}
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
#' @param temperature Optional list with components `T_obs` (observation
#'   temperature in Kelvin), `T_ref` (reference temperature in K), and
#'   `T_A` (Arrhenius temperature in K).  If provided, the correction
#'   factor \eqn{c_T = \exp(T_A/T_{\mathrm{ref}} - T_A/T_{\mathrm{obs}})}
#'   (see [arrhenius()]) is applied **inside the Stan ODE** to all rate
#'   parameters: \eqn{\{p_{Am}\}}, \eqn{[p_M]}, \eqn{v}, \eqn{k_J}
#'   (growth_repro), \eqn{k_d} (debtox).  Parameters that are **not**
#'   rates (\eqn{\kappa}, \eqn{[E_G]}, \eqn{z_w}, \eqn{b_w}) are not
#'   corrected.  The posterior gives estimates at \eqn{T_{\mathrm{ref}}},
#'   making results comparable across experiments at different
#'   temperatures.  The correction is global (single temperature for the
#'   entire experiment).  If `NULL` (default), no correction is applied
#'   (\eqn{c_T = 1}).  The field was renamed from `T` to `T_obs` in
#'   version 0.1.4 to avoid shadowing R's built-in `T` symbol.
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
		cli::cli_abort(c(
			"x" = "{.val individual} model requires exactly 1 individual but data contain {data$n_ind}.",
			"i" = "Filter to one individual, or use {.val hierarchical} for multi-individual analysis."
		))
	}
	if (type == "growth_repro" && !("reproduction" %in% data$endpoints)) {
		cli::cli_abort("growth_repro model requires reproduction data in {.fun bdeb_data}.")
	}
	if (type == "growth_repro" && data$n_ind > 1) {
		cli::cli_abort(c(
			"x" = "{.val growth_repro} model is a single-individual model but data contain {data$n_ind} individuals.",
			"i" = "This model fits one latent DEB process (E, V, R) to one organism.",
			"i" = "Passing multiple individuals would collapse them into one pseudo-organism, producing invalid inference.",
			"i" = "Filter to one individual, or use {.val hierarchical} for multi-individual analysis."
		))
	}
	if (type == "debtox" && is.null(data$concentration)) {
		cli::cli_abort("debtox model requires {.arg concentration} in {.fun bdeb_data}.")
	}

	# Fill defaults
	defaults <- prior_default(type)
	for (nm in names(defaults)) {
		if (is.null(priors[[nm]])) priors[[nm]] <- defaults[[nm]]
	}

	# Validate prior objects
	for (nm in names(priors)) {
		if (!inherits(priors[[nm]], "bdeb_prior")) {
			cli::cli_abort("Prior {.field {nm}} must be a {.cls bdeb_prior} object (use {.fn prior_lognormal}, {.fn prior_beta}, etc.).")
		}
	}

	# Observation models — defaults
	if (is.null(observation$growth))       observation$growth <- obs_normal()
	if (is.null(observation$reproduction)) observation$reproduction <- obs_negbinom()

	# Validate observation objects
	if (!inherits(observation$growth, "bdeb_obs")) {
		cli::cli_abort("{.field observation$growth} must be a {.cls bdeb_obs} object (use {.fn obs_normal}, etc.).")
	}
	if (!inherits(observation$reproduction, "bdeb_obs")) {
		cli::cli_abort("{.field observation$reproduction} must be a {.cls bdeb_obs} object (use {.fn obs_negbinom}, etc.).")
	}

	# Stan model selection
	stan_model_name <- switch(type,
		individual   = "bdeb_individual_growth",
		growth_repro = "bdeb_growth_repro",
		hierarchical = "bdeb_hierarchical_growth",
		debtox       = "bdeb_debtox"
	)

	# Validate temperature if provided
	if (!is.null(temperature)) {
		required_t <- c("T_obs", "T_ref", "T_A")
		missing_t <- setdiff(required_t, names(temperature))
		if (length(missing_t) > 0) {
			cli::cli_abort("Temperature list missing: {.field {missing_t}}. Required: T_obs, T_ref, T_A.")
		}
		for (tn in required_t) {
			val <- temperature[[tn]]
			if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0) {
				cli::cli_abort("Temperature field {.field {tn}} must be a positive finite scalar (Kelvin).")
			}
		}
	}

	# Validate observation model vs type
	g_fam <- observation$growth$family
	if (!(g_fam %in% c("normal", "lognormal", "student_t"))) {
		cli::cli_abort("Growth observation model must be {.val normal}, {.val lognormal}, or {.val student_t}, got {.val {g_fam}}.")
	}
	r_fam <- observation$reproduction$family
	if (!(r_fam %in% c("negbinom", "poisson"))) {
		cli::cli_abort("Reproduction observation model must be {.val negbinom} or {.val poisson}, got {.val {r_fam}}.")
	}

	# Build Stan data
	stan_data <- switch(type,
		individual   = build_stan_data_individual(data, priors, temperature, observation),
		growth_repro = build_stan_data_growth_repro(data, priors, temperature, observation),
		hierarchical = build_stan_data_hierarchical(data, priors, temperature, observation),
		debtox       = build_stan_data_debtox(data, priors, temperature, observation)
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
		cli::cli_alert("  T_obs = {x$temperature$T_obs} K, T_ref = {x$temperature$T_ref} K, T_A = {x$temperature$T_A} K")
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
	if (!is.numeric(nu) || length(nu) != 1 || nu < 1)
		cli::cli_abort("{.arg nu} must be >= 1 (degrees of freedom).")
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

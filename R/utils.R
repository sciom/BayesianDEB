#' Arrhenius Temperature Correction
#'
#' Computes the temperature correction factor for DEB rate parameters
#' based on the Arrhenius relationship (Kooijman, 2010, Eq. 1.2).  In
#' DEB theory all rate parameters (e.g., \eqn{\{p_{Am}\}}, \eqn{[p_M]},
#' \eqn{v}) scale with temperature by the same factor:
#'
#' \deqn{c_T = \exp\!\left(\frac{T_A}{T_{\mathrm{ref}}} - \frac{T_A}{T}\right)}
#'
#' where \eqn{T} and \eqn{T_{\mathrm{ref}}} are in Kelvin and \eqn{T_A}
#' is the Arrhenius temperature (a species-specific constant, typically
#' 6000--12000 K for ectotherms; Kooijman, 2010, Table 8.1).  At
#' \eqn{T = T_{\mathrm{ref}}}, the factor is exactly 1.
#'
#' @param temp Body (or ambient) temperature in Kelvin.  Renamed from
#'   `T` in version 0.1.4 to avoid shadowing R's built-in `T` symbol
#'   (= `TRUE`); pass positionally or use `temp = ...` explicitly.
#' @param T_ref Reference temperature in Kelvin (default 293.15 K = 20 °C).
#' @param T_A Arrhenius temperature in Kelvin (default 8000 K).
#' @return Numeric correction factor (dimensionless, > 0).
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press, Eq. 1.2.
#' \doi{10.1017/CBO9780511805400}
#'
#' @export
#' @examples
#' # Correction at 25 C relative to 20 C reference
#' arrhenius(298.15, T_ref = 293.15, T_A = 8000)  # ~ 1.74
#'
#' # No correction at reference temperature
#' arrhenius(293.15)  # exactly 1
arrhenius <- function(temp, T_ref = 293.15, T_A = 8000) {
	if (!is.numeric(temp) || any(temp <= 0))
		cli::cli_abort("{.arg temp} must be positive (Kelvin).")
	if (!is.numeric(T_ref) || length(T_ref) != 1 || T_ref <= 0)
		cli::cli_abort("{.arg T_ref} must be a positive scalar (Kelvin).")
	if (!is.numeric(T_A) || length(T_A) != 1 || T_A < 0)
		cli::cli_abort("{.arg T_A} must be a non-negative scalar (Kelvin).")
	exp(T_A / T_ref - T_A / temp)
}

#' Compute DEB Energy Fluxes
#'
#' Given current state \eqn{(E, V)} and the core DEB parameters, computes
#' all standard energy fluxes defined by the \eqn{\kappa}-rule
#' (Kooijman, 2010, Eqs. 2.3--2.12):
#'
#' \describe{
#'   \item{\eqn{\dot{p}_A}}{Assimilation: \eqn{f \{p_{Am}\} L^2}.}
#'   \item{\eqn{\dot{p}_C}}{Mobilisation: \eqn{E v L / (E + [E_G] V)}.}
#'   \item{\eqn{\dot{p}_M}}{Somatic maintenance: \eqn{[p_M] V}.}
#'   \item{\eqn{\dot{p}_G}}{Growth: \eqn{\max(\kappa \dot{p}_C - \dot{p}_M, 0)}.}
#'   \item{\eqn{\dot{p}_J}}{Maturity maintenance: \eqn{k_J E_H^p}.}
#'   \item{\eqn{\dot{p}_R}}{Reproduction: \eqn{\max((1-\kappa)\dot{p}_C - \dot{p}_J, 0)}.}
#' }
#'
#' @param E Reserve energy (J).
#' @param V Structural volume (cm\eqn{^3}).
#' @param f Scaled functional response \eqn{f \in [0, 1]}.
#' @param p_Am Surface-area-specific maximum assimilation rate
#'   \eqn{\{p_{Am}\}} (J d\eqn{^{-1}} cm\eqn{^{-2}}).
#' @param p_M Volume-specific somatic maintenance rate \eqn{[p_M]}
#'   (J d\eqn{^{-1}} cm\eqn{^{-3}}).
#' @param kappa Allocation fraction to soma \eqn{\kappa \in (0, 1)}.
#' @param v Energy conductance (cm d\eqn{^{-1}}).
#' @param E_G Specific cost of structure \eqn{[E_G]} (J cm\eqn{^{-3}}).
#' @param k_J Maturity maintenance rate coefficient \eqn{k_J}
#'   (d\eqn{^{-1}}).  Default 0.
#' @param E_Hp Maturity at puberty \eqn{E_H^p} (J).  Default 0.
#' @return Named list with fluxes `p_A`, `p_C`, `p_M`, `p_G`, `p_J`,
#'   `p_R`, structural length `L` (\eqn{V^{1/3}}), and scaled reserve
#'   density `e` (\eqn{E / ([E_m] V + 10^{-12})}).  The \eqn{10^{-12}}
#'   stabilisation prevents division by zero when \eqn{V = 0}; in that
#'   edge case `e` returns a small finite number rather than `Inf` or
#'   `NA`.
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press, Ch. 2.
#' \doi{10.1017/CBO9780511805400}
#' @export
#' @examples
#' # Energy fluxes for a 1 cm^3 organism with typical Eisenia parameters
#' deb_fluxes(E = 1, V = 1, f = 1, p_Am = 5, p_M = 0.5,
#'   kappa = 0.75, v = 0.2, E_G = 4400)
deb_fluxes <- function(E, V, f, p_Am, p_M, kappa, v, E_G,
                       k_J = 0, E_Hp = 0) {
	for (nm in c("E", "V", "f", "p_Am", "p_M", "kappa", "v", "E_G",
	             "k_J", "E_Hp")) {
		x <- get(nm)
		if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
			cli::cli_abort("{.arg {nm}} must be a finite numeric scalar.")
		}
	}
	if (E < 0 || V < 0 || p_Am < 0 || p_M < 0 || v < 0 || E_G < 0 ||
	    k_J < 0 || E_Hp < 0) {
		cli::cli_abort("Energy fluxes require non-negative arguments.")
	}
	if (f < 0 || f > 1) {
		cli::cli_abort("{.arg f} must lie in [0, 1].")
	}
	if (kappa <= 0 || kappa >= 1) {
		cli::cli_abort("{.arg kappa} must lie in (0, 1).")
	}
	L <- V^(1 / 3)
	E_m <- p_Am / v  # maximum reserve density [E_m] = {p_Am}/v

	# Energy fluxes (all in J/d)
	flux_p_A <- f * p_Am * L^2                       # assimilation
	# Mobilisation: p_C = E * v * L / (E + [E_G] * V)
	# Equivalent to [E] * v * [E_G] * L^2 / ([E] + [E_G]) where [E] = E/V.
	# This formulation uses total reserve E (not density [E]) and avoids
	# the need for [E_m].  The 1e-12 stabilises near V = E = 0.
	flux_p_C <- E * v * L / (E + E_G * V + 1e-12)
	flux_p_M <- p_M * V                              # somatic maintenance
	flux_p_J <- k_J * E_Hp                           # maturity maintenance
	flux_p_G <- max(kappa * flux_p_C - flux_p_M, 0)  # growth (energy to structure)
	flux_p_R <- max((1 - kappa) * flux_p_C - flux_p_J, 0)  # reproduction

	list(
		p_A = flux_p_A,
		p_C = flux_p_C,
		p_M = flux_p_M,
		p_G = flux_p_G,
		p_J = flux_p_J,
		p_R = flux_p_R,
		L   = L,
		e   = E / (E_m * V + 1e-12)
	)
}

#' Check that cmdstanr is available
#'
#' Since \pkg{cmdstanr} is listed under Suggests (it is not on CRAN),
#' every function that needs it must call this guard first.
#'
#' @return `TRUE` invisibly if cmdstanr is available; otherwise
#'   throws an informative error.
#' @keywords internal
check_cmdstanr <- function() {
	if (!requireNamespace("cmdstanr", quietly = TRUE)) {
		cli::cli_abort(c(
			"x" = "{.pkg cmdstanr} is required but not installed.",
			"i" = 'Install with: install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))',
			"i" = "Then run: cmdstanr::install_cmdstan()"
		))
	}

	# cmdstanr is installed — now check that CmdStan itself is present
	path <- tryCatch(
		cmdstanr::cmdstan_path(),
		error = function(e) NULL
	)
	if (is.null(path) || !dir.exists(path)) {
		cli::cli_abort(c(
			"x" = "{.pkg cmdstanr} is installed, but CmdStan was not found.",
			"i" = "Install CmdStan with: cmdstanr::install_cmdstan()",
			"i" = "Or set the path manually: cmdstanr::set_cmdstan_path(\"/path/to/cmdstan\")"
		))
	}

	invisible(TRUE)
}

#' Reproducibility Report
#'
#' Prints a concise summary of the software environment used for a
#' BayesianDEB analysis.  Useful for supplementary materials in
#' publications or for diagnosing cross-machine discrepancies.
#'
#' @param fit Optional [bdeb_fit()] object.  If provided, includes
#'   model-specific details (type, chains, iterations, adapt_delta, seed).
#' @return Invisibly returns a named list with all reported information.
#' @export
#' @examples
#' bdeb_session_info()
bdeb_session_info <- function(fit = NULL) {
	info <- list()

	# R version
	info$R_version <- paste0(R.version$major, ".", R.version$minor)
	info$platform  <- R.version$platform

	# Package version
	info$BayesianDEB <- as.character(utils::packageVersion("BayesianDEB"))

	# cmdstanr / CmdStan
	if (requireNamespace("cmdstanr", quietly = TRUE)) {
		info$cmdstanr <- as.character(utils::packageVersion("cmdstanr"))
		cs_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
		if (!is.null(cs_path)) {
			info$cmdstan_version <- tryCatch(
				cmdstanr::cmdstan_version(),
				error = function(e) "unknown"
			)
			info$cmdstan_path <- cs_path
		} else {
			info$cmdstan_version <- "not installed"
		}
	} else {
		info$cmdstanr <- "not installed"
		info$cmdstan_version <- "not installed"
	}

	# Key dependencies
	for (pkg in c("posterior", "loo", "ggplot2", "bayesplot")) {
		info[[pkg]] <- tryCatch(
			as.character(utils::packageVersion(pkg)),
			error = function(e) "not installed"
		)
	}

	# Fit-specific info
	if (!is.null(fit) && inherits(fit, "bdeb_fit")) {
		info$model_type     <- fit$model$type
		info$stan_model     <- fit$model$stan_model_name
		info$chains         <- fit$chains
		info$iter_warmup    <- fit$iter_warmup
		info$iter_sampling  <- fit$iter_sampling
		info$adapt_delta    <- fit$adapt_delta

		# Stan model hash (content fingerprint)
		stan_path <- tryCatch(stan_file(fit$model$stan_model_name),
		                      error = function(e) NULL)
		if (!is.null(stan_path) && file.exists(stan_path) &&
		    requireNamespace("digest", quietly = TRUE)) {
			info$stan_model_hash <- substr(
				digest::digest(readLines(stan_path), algo = "md5"),
				1, 12)
		}
	}

	# Print
	cli::cli_h2("BayesianDEB Reproducibility Report")
	cli::cli_alert_info("R: {info$R_version} ({info$platform})")
	cli::cli_alert_info("BayesianDEB: {info$BayesianDEB}")
	cli::cli_alert_info("cmdstanr: {info$cmdstanr}")
	cli::cli_alert_info("CmdStan: {info$cmdstan_version}")
	cli::cli_alert_info("posterior: {info$posterior}, loo: {info$loo}")

	if (!is.null(fit) && inherits(fit, "bdeb_fit")) {
		cli::cli_h3("Fit Configuration")
		cli::cli_alert("Model: {info$model_type} ({info$stan_model})")
		cli::cli_alert("Chains: {info$chains}, Warmup: {info$iter_warmup}, Sampling: {info$iter_sampling}")
		cli::cli_alert("adapt_delta: {info$adapt_delta}")
		if (!is.null(info$stan_model_hash)) {
			cli::cli_alert("Stan model hash: {info$stan_model_hash}")
		}
	}

	invisible(info)
}

#' Validate finite numeric scalar
#' @keywords internal
assert_finite_scalar <- function(x, name) {
	if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
		cli::cli_abort("{.arg {name}} must be a finite numeric scalar, got {.val {x}}.")
	}
}

#' Validate positive numeric scalar
#' @keywords internal
assert_positive <- function(x, name) {
	if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) {
		cli::cli_abort("{.arg {name}} must be a positive scalar, got {.val {x}}.")
	}
}

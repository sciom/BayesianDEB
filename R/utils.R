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
#' @param T Body (or ambient) temperature in Kelvin.
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
arrhenius <- function(T, T_ref = 293.15, T_A = 8000) {
	exp(T_A / T_ref - T_A / T)
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
#'   density `e` (\eqn{E / ([E_m] V)}).
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press, Ch. 2.
#' \doi{10.1017/CBO9780511805400}
#' @export
deb_fluxes <- function(E, V, f, p_Am, p_M, kappa, v, E_G,
                       k_J = 0, E_Hp = 0) {
	L <- V^(1 / 3)
	E_m <- p_Am / v  # maximum reserve density

	# Fluxes
	flux_p_A <- f * p_Am * L^2
	flux_p_C <- E * v * L / (E + E_G * V + 1e-12)
	flux_p_M <- p_M * V
	flux_p_J <- k_J * E_Hp
	flux_p_G <- max(kappa * flux_p_C - flux_p_M, 0) / E_G * E_G
	flux_p_R <- max((1 - kappa) * flux_p_C - flux_p_J, 0)

	list(
		p_A = flux_p_A,
		p_C = flux_p_C,
		p_M = flux_p_M,
		p_G = max(kappa * flux_p_C - flux_p_M, 0),
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
	invisible(TRUE)
}

#' Validate positive numeric scalar
#' @keywords internal
assert_positive <- function(x, name) {
	if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) {
		cli::cli_abort("{.arg {name}} must be a positive scalar, got {.val {x}}.")
	}
}

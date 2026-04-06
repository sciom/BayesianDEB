#' Simulate DEB Growth Trajectory
#'
#' Forward-integrates the standard 2-state DEB ODE (reserve \eqn{E},
#' structure \eqn{V}) using Euler's method.  This is a standalone
#' simulator independent of Stan — useful for exploring parameter
#' space, generating synthetic data, teaching, and prior predictive
#' checks.
#'
#' @section Numerical note:
#' This uses a fixed-step Euler integrator, which is an approximation.
#' The Stan models use the adaptive BDF solver for fitting.  For
#' visualisation and exploration the Euler integrator is adequate;
#' reduce `dt` for higher accuracy.
#'
#' @param t_max End time (days).
#' @param p_Am Surface-area-specific assimilation rate \eqn{\{p_{Am}\}}
#'   (J d\eqn{^{-1}} cm\eqn{^{-2}}).
#' @param p_M Volume-specific somatic maintenance \eqn{[p_M]}
#'   (J d\eqn{^{-1}} cm\eqn{^{-3}}).
#' @param kappa Allocation fraction to soma.
#' @param v Energy conductance (cm d\eqn{^{-1}}).
#' @param E_G Specific cost of structure (J cm\eqn{^{-3}}).
#' @param E0 Initial reserve density (J cm\eqn{^{-3}}).
#' @param L0 Initial structural length (cm).
#' @param f Scaled functional response \eqn{f \in [0,1]}.  Default 1.
#' @param dt Integration step size (days).  Default 0.5.
#' @return Data frame with columns `time`, `E` (reserve), `V` (volume),
#'   `L` (structural length).
#' @export
#' @examples
#' # Simulate E. fetida growth for 84 days
#' traj <- deb_simulate(t_max = 84, p_Am = 5, p_M = 0.5,
#'   kappa = 0.75, v = 0.2, E_G = 400, E0 = 1, L0 = 0.1)
#' plot(traj$time, traj$L, type = "l", xlab = "Days", ylab = "L (cm)")
deb_simulate <- function(t_max, p_Am, p_M, kappa, v, E_G, E0, L0,
                         f = 1.0, dt = 0.5) {
	assert_positive(t_max, "t_max")
	assert_positive(p_Am, "p_Am")
	assert_positive(p_M, "p_M")
	assert_positive(v, "v")
	assert_positive(E_G, "E_G")
	assert_positive(E0, "E0")
	assert_positive(L0, "L0")

	n <- ceiling(t_max / dt)
	t <- seq(0, t_max, length.out = n + 1)
	E <- V <- numeric(n + 1)
	V[1] <- L0^3
	E[1] <- E0 * V[1]

	for (i in seq_len(n)) {
		L <- V[i]^(1/3)
		pA <- f * p_Am * L^2
		pC <- E[i] * v * L / (E[i] + E_G * V[i] + 1e-12)
		pM <- p_M * V[i]
		dE <- pA - pC
		dV <- (kappa * pC - pM) / E_G
		E[i + 1] <- max(E[i] + dE * dt, 1e-12)
		V[i + 1] <- max(V[i] + dV * dt, 1e-12)
	}
	data.frame(time = t, E = E, V = V, L = V^(1/3))
}

#' Simulate DEBtox Growth Under Toxicant Exposure
#'
#' Forward-integrates the 4-state DEBtox ODE (\eqn{E}, \eqn{V},
#' \eqn{E_R}, \eqn{D_w}) with stress on assimilation.  Standalone
#' simulator independent of Stan.
#'
#' @inheritParams deb_simulate
#' @param k_d Damage recovery rate (d\eqn{^{-1}}).
#' @param z_w No-effect concentration (NEC).
#' @param b_w Effect intensity.
#' @param C_w External toxicant concentration.
#' @return Data frame with columns `time`, `E`, `V`, `L`, `R`
#'   (reproduction buffer), `Dw` (scaled damage).
#' @export
#' @examples
#' traj <- debtox_simulate(t_max = 42, p_Am = 5, p_M = 0.5,
#'   kappa = 0.75, v = 0.2, E_G = 400, E0 = 1, L0 = 0.1,
#'   k_d = 0.3, z_w = 15, b_w = 0.003, C_w = 80)
#' plot(traj$time, traj$L, type = "l")
debtox_simulate <- function(t_max, p_Am, p_M, kappa, v, E_G, E0, L0,
                            k_d, z_w, b_w, C_w, f = 1.0, dt = 0.5) {
	assert_positive(t_max, "t_max")
	assert_positive(p_Am, "p_Am")
	assert_positive(p_M, "p_M")
	assert_positive(v, "v")
	assert_positive(E_G, "E_G")
	assert_positive(E0, "E0")
	assert_positive(L0, "L0")
	assert_positive(k_d, "k_d")

	n <- ceiling(t_max / dt)
	t <- seq(0, t_max, length.out = n + 1)
	E <- V <- Dw <- R <- numeric(n + 1)
	V[1] <- L0^3
	E[1] <- E0 * V[1]

	for (i in seq_len(n)) {
		L <- V[i]^(1/3)
		dDw <- k_d * (max(C_w - z_w, 0) - Dw[i])
		s <- b_w * max(Dw[i], 0)
		pA <- f * p_Am * L^2 * max(1 - s, 0)
		pC <- E[i] * v * L / (E[i] + E_G * V[i] + 1e-12)
		pM <- p_M * V[i]
		dE <- pA - pC
		dV <- (kappa * pC - pM) / E_G
		dR <- max((1 - kappa) * pC, 0)
		E[i + 1]  <- max(E[i] + dE * dt, 1e-12)
		V[i + 1]  <- max(V[i] + dV * dt, 1e-12)
		Dw[i + 1] <- max(Dw[i] + dDw * dt, 0)
		R[i + 1]  <- R[i] + dR * dt
	}
	data.frame(time = t, E = E, V = V, L = V^(1/3), R = R, Dw = Dw)
}

#' @keywords internal
"_PACKAGE"

#' BayesianDEB: Bayesian Dynamic Energy Budget Modelling
#'
#' Embeds the standard Dynamic Energy Budget (DEB) model of Kooijman (2010)
#' within a Bayesian state-space framework, using Stan (Carpenter et al., 2017)
#' for Hamiltonian Monte Carlo inference.  The DEB ordinary differential
#' equation system describes how organisms allocate assimilated energy to
#' maintenance, growth, and reproduction according to the \eqn{\kappa}-rule.
#' This package wraps pre-compiled Stan programs with a declarative R
#' interface so that users can fit DEB models without writing Stan code.
#'
#' @section Models:
#' Four model types are implemented, covering the most common DEB
#' applications in ecology and ecotoxicology:
#' \describe{
#'   \item{`"individual"`}{Two-state model (reserve \eqn{E}, structure
#'     \eqn{V}) for a single organism.  The ODE follows Eqs. 2.4 and 2.6
#'     of Kooijman (2010, Ch. 2).}
#'   \item{`"growth_repro"`}{Three-state model (\eqn{E}, \eqn{V},
#'     reproduction buffer \eqn{E_R}) with negative-binomial observation
#'     model for offspring counts.}
#'   \item{`"hierarchical"`}{Two-state model with lognormal random effects
#'     on the surface-area-specific assimilation rate \eqn{\{p_{Am}\}},
#'     using the non-centred parameterisation of Betancourt & Girolami
#'     (2015) to avoid pathological funnel geometry.}
#'   \item{`"debtox"`}{Four-state TKTD model (\eqn{E}, \eqn{V},
#'     \eqn{E_R}, scaled damage \eqn{D_w}) following the DEBtox framework
#'     of Jager et al. (2006).  Stress on assimilation is currently
#'     implemented; the damage dynamics follow
#'     \eqn{dD_w/dt = k_d(\max(C_w - z_w, 0) - D_w)}.}
#' }
#'
#' @section Workflow:
#' The recommended workflow follows the iterative Bayesian modelling cycle
#' described in Gelman et al. (2013, Ch. 6):
#' \enumerate{
#'   \item Prepare data with [bdeb_data()].
#'   \item Specify model and priors with [bdeb_model()].
#'   \item Fit via MCMC with [bdeb_fit()].
#'   \item Check convergence with [bdeb_diagnose()]
#'     (\eqn{\hat{R}}, ESS, divergences; Vehtari et al., 2021).
#'   \item Perform posterior predictive checks with [bdeb_ppc()]
#'     (Gelman et al., 2013, Ch. 6).
#'   \item Compute derived quantities with [bdeb_derived()].
#'   \item Compare models via [bdeb_loo()] (individual and
#'     growth_repro models only).
#'   \item Iterate: revise priors or model structure as needed.
#' }
#'
#' @section Key DEB equations:
#' \deqn{\dot{p}_A = f \{p_{Am}\} L^2}
#' \deqn{\dot{p}_C = \frac{E \cdot v \cdot L}{E + [E_G] V}}
#' \deqn{dE/dt = \dot{p}_A - \dot{p}_C}
#' \deqn{dV/dt = (\kappa \dot{p}_C - [p_M] V) / [E_G]}
#' where \eqn{L = V^{1/3}} is structural length, \eqn{f} is the scaled
#' functional response, and all symbols follow the DEB notation of
#' Kooijman (2010, Table 1.1).
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press.
#' \doi{10.1017/CBO9780511805400}
#'
#' Carpenter, B., Gelman, A., Hoffman, M.D., Lee, D., Goodrich, B.,
#' Betancourt, M., Brubaker, M.A., Guo, J., Li, P. and Riddell, A.
#' (2017). Stan: A probabilistic programming language. *Journal of
#' Statistical Software*, 76(1), 1--32. \doi{10.18637/jss.v076.i01}
#'
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A.
#' and Rubin, D.B. (2013). *Bayesian Data Analysis*. 3rd edition.
#' Chapman & Hall/CRC.
#'
#' Jager, T., Heugens, E.H.W. and Kooijman, S.A.L.M. (2006). Making
#' sense of ecotoxicological test results: towards application of
#' process-based models. *Ecotoxicology*, 15(3), 305--314.
#' \doi{10.1007/s10646-006-0060-x}
#'
#' Betancourt, M. and Girolami, M. (2015). Hamiltonian Monte Carlo for
#' hierarchical models. In: Upadhyay, S.K. et al. (eds) *Current Trends
#' in Bayesian Methodology with Applications*. CRC Press, pp. 79--101.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Bürkner, P.-C. (2021). Rank-normalization, folding, and localization:
#' an improved \eqn{\hat{R}} for assessing convergence of MCMC.
#' *Bayesian Analysis*, 16(2), 667--718. \doi{10.1214/20-BA1221}
#'
#' @section Numerical layers:
#' The package uses two distinct numerical engines:
#' \describe{
#'   \item{Stan inference (exact)}{`bdeb_fit()`, `bdeb_diagnose()`,
#'     `bdeb_summary()`, `bdeb_ec50()`, `bdeb_loo()`, `bdeb_ppc()`
#'     for individual/growth_repro.  These use Stan's BDF stiff ODE
#'     solver with adaptive step size and tolerances \eqn{10^{-6}}.
#'     Posteriors, derived quantities, and log-likelihoods from this
#'     layer are publication-grade.}
#'   \item{R-side simulation}{`bdeb_prior_predictive()`,
#'     `bdeb_predict(newdata=...)`,
#'     `plot(fit, type="trajectory")` for hierarchical/debtox,
#'     `plot_dose_response()`.  These use the LSODA solver from
#'     \pkg{deSolve} (adaptive step size, tolerances \eqn{10^{-6}}),
#'     matching Stan's BDF solver.  They are **visualisation and
#'     exploration tools**, not exact inference.  For quantitative
#'     results, use the Stan-based functions.}
#' }
#'
#' @section Lifecycle:
#' \describe{
#'   \item{Stable}{`bdeb_data()`, `bdeb_model()`, `bdeb_fit()`,
#'     `bdeb_diagnose()`, `bdeb_summary()`, `bdeb_derived()`,
#'     `bdeb_ppc()`, `bdeb_predict()`, `bdeb_loo()`, `bdeb_ec50()`,
#'     `bdeb_tox()`, `bdeb_prior_predictive()`, `plot_dose_response()`,
#'     `coef()`, all prior and observation model constructors,
#'     `arrhenius()`, `deb_fluxes()`, `repro_to_intervals()`,
#'     `bdeb_session_info()`.}
#'   \item{Stable models}{`"individual"`, `"growth_repro"`,
#'     `"hierarchical"`.}
#'   \item{Beta models}{`"debtox"` (group-level, 1 ODE per concentration;
#'     hierarchical individual-level DEBtox is planned).}
#'   \item{Planned}{Survival endpoint, per-observation temperature,
#'     DEBtox stress on maintenance/growth cost, hierarchical DEBtox
#'     (individual-level TKTD), weight endpoint with shape coefficient
#'     estimation, wider maturity dynamics, real-data benchmark dataset.}
#' }
#'
#' @name BayesianDEB-package
#' @aliases BayesianDEB
NULL

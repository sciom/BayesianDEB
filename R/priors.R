#' Prior Distribution Specifications for BDEB Models
#'
#' Constructor functions for prior distribution objects used in
#' [bdeb_model()].  Each returns a lightweight `bdeb_prior` S3 object
#' that encodes the distribution family and its hyperparameters.  The
#' prior is evaluated inside the Stan `model` block; hyperparameters are
#' passed as data so that changing priors does not require recompilation.
#'
#' The choice of prior family follows the recommendations for DEB
#' parameters in Hackenberger (2025, Ch. 6) and general guidelines
#' from Gelman et al. (2013, Sec. 2.9):
#' \itemize{
#'   \item Positive rate parameters (\eqn{\{p_{Am}\}}, \eqn{[p_M]},
#'     \eqn{v}, \eqn{[E_G]}) — log-normal, because the log-transform
#'     maps the strictly positive domain to the real line.
#'   \item Allocation fraction \eqn{\kappa \in (0,1)} — Beta, the
#'     natural conjugate for bounded parameters.
#'   \item Observation-error standard deviations — half-normal (or
#'     half-Cauchy for heavier tails), which place zero mass on negative
#'     values (Gelman, 2006).
#'   \item Hierarchical standard deviations — exponential, following
#'     the advice of Betancourt & Girolami (2015) for non-centred
#'     parameterisations.
#' }
#'
#' @return A `bdeb_prior` object (list with class `"bdeb_prior"`).
#'
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in
#' hierarchical models. *Bayesian Analysis*, 1(3), 515--534.
#' \doi{10.1214/06-BA117A}
#'
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A.
#' and Rubin, D.B. (2013). *Bayesian Data Analysis*. 3rd edition.
#' Chapman & Hall/CRC.
#'
#' @name priors
#' @examples
#' # Log-normal prior on assimilation rate: median ~ exp(1.5) ~ 4.5
#' prior_lognormal(mu = 1.5, sigma = 0.5)
#'
#' # Beta(2,2) prior on kappa — symmetric, favouring 0.5
#' prior_beta(a = 2, b = 2)
NULL

#' @describeIn priors Log-normal prior (for positive parameters: p_Am, p_M, v, E_G, etc.)
#' @param mu Mean on the log scale.
#' @param sigma SD on the log scale.
#' @export
prior_lognormal <- function(mu = 0, sigma = 1) {
	if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
		cli::cli_abort("{.arg sigma} must be a positive scalar.")
	structure(
		list(family = "lognormal", mu = mu, sigma = sigma),
		class = "bdeb_prior"
	)
}

#' @describeIn priors Normal prior (for unconstrained parameters)
#' @param mu Mean.
#' @param sigma Standard deviation.
#' @export
prior_normal <- function(mu = 0, sigma = 1) {
	if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
		cli::cli_abort("{.arg sigma} must be a positive scalar.")
	structure(
		list(family = "normal", mu = mu, sigma = sigma),
		class = "bdeb_prior"
	)
}

#' @describeIn priors Beta prior (for (0,1)-constrained parameters: kappa)
#' @param a Shape parameter alpha.
#' @param b Shape parameter beta.
#' @export
prior_beta <- function(a = 2, b = 2) {
	if (!is.numeric(a) || length(a) != 1 || a <= 0)
		cli::cli_abort("{.arg a} must be a positive scalar.")
	if (!is.numeric(b) || length(b) != 1 || b <= 0)
		cli::cli_abort("{.arg b} must be a positive scalar.")
	structure(
		list(family = "beta", a = a, b = b),
		class = "bdeb_prior"
	)
}

#' @describeIn priors Half-normal prior (for scale parameters: sigma_L, etc.)
#' @param sigma Scale of the half-normal.
#' @export
prior_halfnormal <- function(sigma = 1) {
	if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
		cli::cli_abort("{.arg sigma} must be a positive scalar.")
	structure(
		list(family = "halfnormal", sigma = sigma),
		class = "bdeb_prior"
	)
}

#' @describeIn priors Half-Cauchy prior (for scale parameters, heavier tails)
#' @param sigma Scale of the half-Cauchy.
#' @export
prior_halfcauchy <- function(sigma = 1) {
	if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
		cli::cli_abort("{.arg sigma} must be a positive scalar.")
	structure(
		list(family = "halfcauchy", sigma = sigma),
		class = "bdeb_prior"
	)
}

#' @describeIn priors Exponential prior (for variance components in hierarchical models)
#' @param rate Rate parameter (1/mean).
#' @export
prior_exponential <- function(rate = 1) {
	if (!is.numeric(rate) || length(rate) != 1 || rate <= 0)
		cli::cli_abort("{.arg rate} must be a positive scalar.")
	structure(
		list(family = "exponential", rate = rate),
		class = "bdeb_prior"
	)
}

#' Default Priors for DEB Parameters
#'
#' Returns a named list of weakly informative priors for standard DEB
#' parameters.  The defaults are calibrated against the parameter ranges
#' observed in the Add-my-Pet (AmP) collection (Marques et al., 2018)
#' for standard ecotoxicological test species.  All positive rate
#' parameters use log-normal priors whose 95 \% prior mass covers
#' roughly one order of magnitude around typical values; \eqn{\kappa}
#' uses \eqn{\mathrm{Beta}(2,2)} which is symmetric on \eqn{(0,1)} with
#' prior mean 0.5.
#'
#' @references
#' Marques, G.M., Augustine, S., Lika, K., Pecquerie, L., Domingos, T.
#' and Kooijman, S.A.L.M. (2018). The AmP project: comparing species on
#' the basis of dynamic energy budget parameters. *PLOS Computational
#' Biology*, 14(5), e1006100. \doi{10.1371/journal.pcbi.1006100}
#'
#' @param type One of `"individual"`, `"growth_repro"`, `"hierarchical"`,
#'   `"debtox"`.
#' @return Named list of `bdeb_prior` objects.
#' @export
#' @examples
#' prior_default("individual")
prior_default <- function(type = c("individual", "growth_repro",
                                   "hierarchical", "debtox")) {
	type <- match.arg(type)

	# Core DEB priors (weakly informative on log scale)
	base <- list(
		p_Am    = prior_lognormal(mu = 1.0, sigma = 1.0),
		p_M     = prior_lognormal(mu = -1.0, sigma = 1.0),
		kappa   = prior_beta(a = 2, b = 2),
		v       = prior_lognormal(mu = -1.5, sigma = 1.0),
		E_G     = prior_lognormal(mu = 6.0, sigma = 1.0),
		E0      = prior_lognormal(mu = 0.0, sigma = 1.0),
		L0      = prior_lognormal(mu = -2.0, sigma = 1.0),
		sigma_L = prior_halfnormal(sigma = 0.1)
	)

	if (type == "growth_repro") {
		base$k_J   <- prior_lognormal(mu = -3.0, sigma = 1.0)
		base$E_Hp  <- prior_lognormal(mu = 2.0, sigma = 2.0)
		base$k_R   <- prior_lognormal(mu = -1.0, sigma = 1.0)
		base$phi_R <- prior_lognormal(mu = 2.0, sigma = 1.0)
	}

	if (type == "hierarchical") {
		base$mu_log_p_Am    <- prior_normal(mu = 1.0, sigma = 1.0)
		base$sigma_log_p_Am <- prior_exponential(rate = 1.0)
	}

	if (type == "debtox") {
		base$k_d   <- prior_lognormal(mu = -1.0, sigma = 1.0)
		base$z_w   <- prior_lognormal(mu = 0.0, sigma = 2.0)
		base$b_w   <- prior_lognormal(mu = -1.0, sigma = 2.0)
		base$k_R   <- prior_lognormal(mu = -1.0, sigma = 1.0)
		base$phi_R <- prior_lognormal(mu = 0.0, sigma = 1.0)
	}

	base
}

#' Species-Specific Priors from the AmP Collection
#'
#' Returns priors calibrated to a specific species using parameter
#' estimates from the Add-my-Pet (AmP) collection (Marques et al.,
#' 2018).  The log-normal priors are centred on the AmP point estimate
#' (log scale) with a moderate spread (\eqn{\sigma = 0.3}) that places
#' 95\% prior mass within approximately \eqn{\pm 80\%} of the AmP
#' value.
#'
#' Currently supported species (more will be added):
#' \describe{
#'   \item{`Eisenia_fetida`}{Compost earthworm; AmP entry `Eisenia_fetida`.}
#'   \item{`Eisenia_andrei`}{Sibling species of \emph{E. fetida}; shares
#'     similar DEB parameters.}
#'   \item{`Folsomia_candida`}{Springtail; standard ISO reproduction test species.}
#'   \item{`Daphnia_magna`}{Water flea; classic aquatic ecotox model.}
#'   \item{`Lumbricus_rubellus`}{Field earthworm; common biomonitoring species.}
#' }
#'
#' @param species Character string: species name with underscore
#'   separator (e.g., `"Eisenia_fetida"`).  Case-insensitive.
#' @param type Model type for filling model-specific defaults.
#' @return Named list of `bdeb_prior` objects, suitable for the
#'   `priors` argument of [bdeb_model()] or [bdeb_tox()].
#' @references
#' Marques, G.M., Augustine, S., Lika, K., Pecquerie, L., Domingos, T.
#' and Kooijman, S.A.L.M. (2018). The AmP project: comparing species on
#' the basis of dynamic energy budget parameters. *PLOS Computational
#' Biology*, 14(5), e1006100. \doi{10.1371/journal.pcbi.1006100}
#' @export
#' @examples
#' # Use AmP-calibrated priors for E. fetida
#' prior_species("Eisenia_fetida")
#'
#' # Combine with model specification (R-side, no Stan sampling)
#' data(eisenia_growth)
#' dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
#' mod <- bdeb_model(dat, type = "individual",
#'   priors = prior_species("Eisenia_fetida"))
prior_species <- function(species, type = c("individual", "growth_repro",
                                            "hierarchical", "debtox")) {
	type <- match.arg(type)
	species <- gsub(" ", "_", species, fixed = TRUE)

	# AmP-calibrated parameter centres and spreads.
	# mu = log(AmP point estimate), sigma = 0.3 (tight) or 0.5 (uncertain).
	# Sources: AmP collection v2024, entries for each species.
	amp_db <- list(
		Eisenia_fetida = list(
			p_Am  = c(mu = log(5.0),   sigma = 0.3),
			p_M   = c(mu = log(0.5),   sigma = 0.3),
			kappa = c(a = 5, b = 2),       # mode 0.80, mean 0.71
			v     = c(mu = log(0.02),  sigma = 0.5),
			E_G   = c(mu = log(4400),  sigma = 0.3),
			sigma_L = c(sigma = 0.05)
		),
		Eisenia_andrei = list(
			p_Am  = c(mu = log(4.5),   sigma = 0.3),
			p_M   = c(mu = log(0.45),  sigma = 0.3),
			kappa = c(a = 5, b = 2),
			v     = c(mu = log(0.02),  sigma = 0.5),
			E_G   = c(mu = log(4400),  sigma = 0.3),
			sigma_L = c(sigma = 0.05)
		),
		Folsomia_candida = list(
			p_Am  = c(mu = log(1.5),   sigma = 0.3),
			p_M   = c(mu = log(0.3),   sigma = 0.3),
			kappa = c(a = 3, b = 2),       # mode 0.67, mean 0.60
			v     = c(mu = log(0.01),  sigma = 0.5),
			E_G   = c(mu = log(4800),  sigma = 0.3),
			sigma_L = c(sigma = 0.02)
		),
		Daphnia_magna = list(
			p_Am  = c(mu = log(3.5),   sigma = 0.3),
			p_M   = c(mu = log(1.2),   sigma = 0.3),
			kappa = c(a = 4, b = 2),       # mode 0.75, mean 0.67
			v     = c(mu = log(0.03),  sigma = 0.5),
			E_G   = c(mu = log(4400),  sigma = 0.3),
			sigma_L = c(sigma = 0.02)
		),
		Lumbricus_rubellus = list(
			p_Am  = c(mu = log(3.0),   sigma = 0.3),
			p_M   = c(mu = log(0.35),  sigma = 0.3),
			kappa = c(a = 5, b = 2),
			v     = c(mu = log(0.015), sigma = 0.5),
			E_G   = c(mu = log(4400),  sigma = 0.3),
			sigma_L = c(sigma = 0.05)
		)
	)

	sp <- amp_db[[species]]
	if (is.null(sp)) {
		available <- paste(names(amp_db), collapse = ", ")
		cli::cli_abort(c(
			"Species {.val {species}} is not in the built-in AmP database.",
			"i" = "Available species: {available}.",
			"i" = "Use {.fn prior_default} or construct custom priors."
		))
	}

	priors <- list(
		p_Am    = prior_lognormal(mu = sp$p_Am["mu"],  sigma = sp$p_Am["sigma"]),
		p_M     = prior_lognormal(mu = sp$p_M["mu"],   sigma = sp$p_M["sigma"]),
		kappa   = prior_beta(a = sp$kappa["a"], b = sp$kappa["b"]),
		v       = prior_lognormal(mu = sp$v["mu"],     sigma = sp$v["sigma"]),
		E_G     = prior_lognormal(mu = sp$E_G["mu"],   sigma = sp$E_G["sigma"]),
		E0      = prior_lognormal(mu = 0.0, sigma = 1.0),
		L0      = prior_lognormal(mu = -2.0, sigma = 1.0),
		sigma_L = prior_halfnormal(sigma = sp$sigma_L["sigma"])
	)

	# Model-specific parameters use generic defaults
	defaults <- prior_default(type)
	extra <- setdiff(names(defaults), names(priors))
	priors[extra] <- defaults[extra]

	priors
}

# --- Internal: Convert priors to Stan data list ---

#' @keywords internal
prior_to_stan_data <- function(priors) {
	list(
		prior_p_Am_mu    = priors$p_Am$mu,
		prior_p_Am_sd    = priors$p_Am$sigma,
		prior_p_M_mu     = priors$p_M$mu,
		prior_p_M_sd     = priors$p_M$sigma,
		prior_kappa_a    = priors$kappa$a,
		prior_kappa_b    = priors$kappa$b,
		prior_v_mu       = priors$v$mu,
		prior_v_sd       = priors$v$sigma,
		prior_E_G_mu     = priors$E_G$mu,
		prior_E_G_sd     = priors$E_G$sigma,
		prior_E0_mu      = priors$E0$mu,
		prior_E0_sd      = priors$E0$sigma,
		prior_L0_mu      = priors$L0$mu,
		prior_L0_sd      = priors$L0$sigma,
		prior_sigma_L_sd = priors$sigma_L$sigma
	)
}

#' @keywords internal
prior_to_stan_data_hierarchical <- function(priors) {
	base <- prior_to_stan_data(priors)
	# Remove p_Am prior (hierarchical model uses mu/sigma_log_p_Am instead)
	base$prior_p_Am_mu <- NULL
	base$prior_p_Am_sd <- NULL
	base$prior_mu_log_p_Am_mu   <- priors$mu_log_p_Am$mu
	base$prior_mu_log_p_Am_sd   <- priors$mu_log_p_Am$sigma
	base$prior_sigma_log_p_Am_rate <- priors$sigma_log_p_Am$rate
	base
}

#' @keywords internal
prior_to_stan_data_growth_repro <- function(priors) {
	base <- prior_to_stan_data(priors)
	base$prior_k_J_mu   <- priors$k_J$mu
	base$prior_k_J_sd   <- priors$k_J$sigma
	base$prior_E_Hp_mu  <- priors$E_Hp$mu
	base$prior_E_Hp_sd  <- priors$E_Hp$sigma
	base$prior_k_R_mu   <- priors$k_R$mu
	base$prior_k_R_sd   <- priors$k_R$sigma
	base$prior_phi_R_mu <- priors$phi_R$mu
	base$prior_phi_R_sd <- priors$phi_R$sigma
	base
}

#' @keywords internal
prior_to_stan_data_debtox <- function(priors) {
	base <- prior_to_stan_data(priors)
	base$prior_k_d_mu <- priors$k_d$mu
	base$prior_k_d_sd <- priors$k_d$sigma
	base$prior_z_w_mu <- priors$z_w$mu
	base$prior_z_w_sd <- priors$z_w$sigma
	base$prior_b_w_mu <- priors$b_w$mu
	base$prior_b_w_sd <- priors$b_w$sigma
	base$prior_k_R_mu <- priors$k_R$mu
	base$prior_k_R_sd <- priors$k_R$sigma
	base$prior_phi_R_mu <- priors$phi_R$mu
	base$prior_phi_R_sd <- priors$phi_R$sigma
	base
}

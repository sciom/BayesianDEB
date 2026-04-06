#' Fit a BDEB Model via Hamiltonian Monte Carlo
#'
#' Compiles the bundled Stan program (cached after the first call) and
#' runs the No-U-Turn Sampler (NUTS; Hoffman & Gelman, 2014) via
#' \pkg{cmdstanr}.  The Stan ODE system is solved at each leapfrog step
#' using the BDF stiff solver (`ode_bdf`) with absolute and relative
#' tolerances of \eqn{10^{-6}}.
#'
#' **Tuning guidance.** If [bdeb_diagnose()] reports divergent transitions,
#' increase `adapt_delta` toward 1.0 (e.g., 0.95 or 0.99).  This reduces
#' the step size, trading speed for geometric fidelity of the sampler.  If
#' the maximum treedepth is frequently saturated, increase
#' `max_treedepth` (e.g., 12 or 15).  For hierarchical models, starting
#' values can matter; the non-centred parameterisation used in
#' `bdeb_hierarchical_growth.stan` should suffice in most cases.
#'
#' @param model A [bdeb_model()] object.
#' @param chains Number of independent MCMC chains. Default 4 (the minimum
#'   recommended by Vehtari et al., 2021, for reliable \eqn{\hat{R}}).
#' @param iter_warmup Number of warmup (adaptation) iterations per chain.
#'   Default 1000.  Stan uses dual-averaging to tune the step size and
#'   diagonal mass matrix during warmup.
#' @param iter_sampling Number of post-warmup sampling iterations per chain.
#'   Default 1000 (yielding 4000 total draws with 4 chains).
#' @param adapt_delta Target Metropolis acceptance probability for NUTS.
#'   Default 0.8.  Increase toward 1.0 to reduce divergences at the cost
#'   of smaller step sizes and longer runtime.
#' @param max_treedepth Maximum binary-tree depth for NUTS.  Default 10
#'   (i.e., up to \eqn{2^{10} = 1024} leapfrog steps per transition).
#' @param seed Integer random seed for full reproducibility.
#' @param parallel_chains Number of chains to run in parallel.
#'   Default `min(chains, detectCores() - 1)`.
#' @param refresh How often to print sampling progress (iterations).
#'   Default 200.  Set to 0 for silent operation.
#' @param ... Additional arguments forwarded to `CmdStanModel$sample()`.
#' @return A `bdeb_fit` object containing the `CmdStanMCMC` result,
#'   the model specification, and sampling metadata.
#'
#' @references
#' Hoffman, M.D. and Gelman, A. (2014). The No-U-Turn Sampler:
#' adaptively setting path lengths in Hamiltonian Monte Carlo.
#' *Journal of Machine Learning Research*, 15(47), 1593--1623.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#' Bürkner, P.-C. (2021). Rank-normalization, folding, and localization:
#' an improved \eqn{\hat{R}} for assessing convergence of MCMC.
#' *Bayesian Analysis*, 16(2), 667--718. \doi{10.1214/20-BA1221}
#' @export
#' @examples
#' \dontrun{
#' dat <- bdeb_data(growth = growth_df)
#' mod <- bdeb_model(dat, type = "individual")
#' fit <- bdeb_fit(mod, chains = 2, iter_sampling = 500)
#' }
bdeb_fit <- function(model,
                     chains = 4,
                     iter_warmup = 1000,
                     iter_sampling = 1000,
                     adapt_delta = 0.8,
                     max_treedepth = 10,
                     seed = NULL,
                     parallel_chains = NULL,
                     refresh = 200,
                     ...) {

	if (!inherits(model, "bdeb_model")) {
		cli::cli_abort("{.arg model} must be a {.cls bdeb_model} object.")
	}

	check_cmdstanr()

	if (is.null(parallel_chains)) {
		n_cores <- tryCatch(parallel::detectCores(), error = function(e) 2L)
		parallel_chains <- min(chains, max(1L, n_cores - 1L))
	}

	# Compile Stan model
	stan_path <- stan_file(model$stan_model_name)
	cli::cli_alert_info("Compiling Stan model: {.file {model$stan_model_name}}")

	stan_mod <- cmdstanr::cmdstan_model(
		stan_file = stan_path,
		cpp_options = list(stan_threads = FALSE)
	)

	# Sample
	cli::cli_alert_info("Running MCMC ({chains} chains, {iter_sampling} iterations each)")

	fit <- stan_mod$sample(
		data            = model$stan_data,
		chains          = chains,
		parallel_chains = parallel_chains,
		iter_warmup     = iter_warmup,
		iter_sampling   = iter_sampling,
		adapt_delta     = adapt_delta,
		max_treedepth   = max_treedepth,
		seed            = seed,
		refresh         = refresh,
		...
	)

	# Construct result
	out <- list(
		fit         = fit,
		model       = model,
		stan_model  = stan_mod,
		chains      = chains,
		iter_warmup = iter_warmup,
		iter_sampling = iter_sampling,
		adapt_delta = adapt_delta
	)

	structure(out, class = "bdeb_fit")
}

#' @return The input object, invisibly.
#' @export
print.bdeb_fit <- function(x, ...) {
	cli::cli_h2("BDEB Fit")
	cli::cli_alert_info("Model type: {x$model$type}")
	cli::cli_alert_info("Chains: {x$chains}, Warmup: {x$iter_warmup}, Sampling: {x$iter_sampling}")

	# Quick diagnostics
	diag <- x$fit$diagnostic_summary(quiet = TRUE)
	n_div <- sum(diag$num_divergent)
	n_tree <- sum(diag$num_max_treedepth)

	if (n_div > 0) {
		cli::cli_alert_warning("Divergent transitions: {n_div}")
	} else {
		cli::cli_alert_success("No divergent transitions")
	}
	if (n_tree > 0) {
		cli::cli_alert_warning("Max treedepth hit: {n_tree} times")
	}

	invisible(x)
}

#' @return A `posterior::draws_summary` data frame (see [bdeb_summary()]).
#' @export
summary.bdeb_fit <- function(object, ...) {
	bdeb_summary(object, ...)
}

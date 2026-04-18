#' Fit a BDEB Model via Hamiltonian Monte Carlo
#'
#' Compiles the bundled Stan program via \pkg{cmdstanr} (which handles
#' caching internally based on the Stan source hash) and runs the
#' No-U-Turn Sampler (NUTS; Hoffman & Gelman, 2014).  The Stan ODE system is solved at each leapfrog step
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
#' @param threads_per_chain Number of threads per chain for within-chain
#'   parallelism via Stan's `reduce_sum`.  Default `NULL` (no threading;
#'   the model runs sequentially within each chain).  Set to e.g.
#'   4 for the `"hierarchical"` or `"debtox"` models to distribute
#'   per-individual / per-group ODE solves across threads.  When
#'   threading is active, adjust `parallel_chains` so that
#'   `chains * threads_per_chain` does not exceed available cores.
#'   Has no effect on `"individual"` or `"growth_repro"` models
#'   (single ODE solve, nothing to distribute).
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
#' # \dontrun{} because bdeb_fit() requires the external CmdStan toolchain
#' # (not on CRAN) and a single fit takes > 30 seconds (Stan compilation
#' # + MCMC).  Users can run this manually after `cmdstanr::install_cmdstan()`.
#' \dontrun{
#' data(eisenia_growth)
#' dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
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
                     threads_per_chain = NULL,
                     refresh = 200,
                     ...) {

	if (!inherits(model, "bdeb_model")) {
		cli::cli_abort("{.arg model} must be a {.cls bdeb_model} object.")
	}
	if (!is.numeric(chains) || chains < 1)
		cli::cli_abort("{.arg chains} must be >= 1.")
	if (!is.numeric(iter_warmup) || iter_warmup < 0)
		cli::cli_abort("{.arg iter_warmup} must be >= 0.")
	if (!is.numeric(iter_sampling) || iter_sampling < 1)
		cli::cli_abort("{.arg iter_sampling} must be >= 1.")
	if (!is.numeric(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1)
		cli::cli_abort("{.arg adapt_delta} must be in (0, 1).")
	if (!is.numeric(max_treedepth) || max_treedepth < 1)
		cli::cli_abort("{.arg max_treedepth} must be >= 1.")
	if (!is.null(threads_per_chain) &&
	    (!is.numeric(threads_per_chain) || threads_per_chain < 1))
		cli::cli_abort("{.arg threads_per_chain} must be >= 1 or NULL.")

	check_cmdstanr()

	# Determine threading
	use_threads <- !is.null(threads_per_chain) && threads_per_chain > 1L
	if (use_threads && model$type %in% c("individual", "growth_repro")) {
		cli::cli_alert_info(
			"Threading has no effect on {.val {model$type}} models (single ODE solve). Ignoring."
		)
		use_threads <- FALSE
	}

	n_cores <- tryCatch(parallel::detectCores(), error = function(e) 2L)

	if (is.null(parallel_chains)) {
		if (use_threads) {
			# Leave room for threads: total = parallel_chains * threads_per_chain
			parallel_chains <- max(1L, min(chains, n_cores %/% threads_per_chain))
		} else {
			parallel_chains <- min(chains, max(1L, n_cores - 1L))
		}
	}

	# Compile Stan model
	stan_path <- stan_file(model$stan_model_name)
	cpp_opts <- list(stan_threads = use_threads)

	if (use_threads) {
		cli::cli_alert_info(
			"Compiling Stan model with threading: {.file {model$stan_model_name}} ({threads_per_chain} threads/chain)"
		)
	} else {
		cli::cli_alert_info("Compiling Stan model: {.file {model$stan_model_name}}")
	}

	stan_mod <- tryCatch(
		cmdstanr::cmdstan_model(stan_file = stan_path, cpp_options = cpp_opts),
		error = function(e) {
			cli::cli_abort(c(
				"x" = "Stan model compilation failed.",
				"i" = "Check your CmdStan installation and C++ toolchain.",
				"x" = conditionMessage(e)
			))
		}
	)

	# Sample
	if (use_threads) {
		cli::cli_alert_info(
			"Running MCMC ({chains} chains, {iter_sampling} iter, {threads_per_chain} threads/chain)"
		)
	} else {
		cli::cli_alert_info("Running MCMC ({chains} chains, {iter_sampling} iterations each)")
	}

	sample_args <- list(
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

	if (use_threads) {
		sample_args$threads_per_chain <- as.integer(threads_per_chain)
	}

	fit <- tryCatch(
		do.call(stan_mod$sample, sample_args),
		error = function(e) {
			cli::cli_abort(c(
				"x" = "MCMC sampling failed.",
				"i" = "Try increasing {.arg adapt_delta} or checking initial values.",
				"x" = conditionMessage(e)
			))
		}
	)

	# Construct result with reproducibility metadata
	out <- list(
		fit               = fit,
		model             = model,
		stan_model        = stan_mod,
		chains            = as.integer(chains),
		iter_warmup       = as.integer(iter_warmup),
		iter_sampling     = as.integer(iter_sampling),
		adapt_delta       = adapt_delta,
		max_treedepth     = as.integer(max_treedepth),
		seed              = seed,
		threads_per_chain = if (use_threads) as.integer(threads_per_chain) else 1L,
		package_version   = as.character(utils::packageVersion("BayesianDEB")),
		cmdstanr_version  = tryCatch(
			as.character(utils::packageVersion("cmdstanr")),
			error = function(e) NA_character_),
		timestamp         = Sys.time()
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

#' Extract Point Estimates from a BDEB Fit
#'
#' Returns posterior medians of all model parameters as a named numeric
#' vector, consistent with the S3 [coef()] convention.
#'
#' @param object A [bdeb_fit()] object.
#' @param type One of `"median"` (default) or `"mean"`.
#' @param ... Ignored.
#' @return Named numeric vector of point estimates.
#' @export
coef.bdeb_fit <- function(object, type = c("median", "mean"), ...) {
	type <- match.arg(type)
	draws <- posterior::as_draws_df(object$fit$draws())
	all_vars <- posterior::variables(draws)
	pars <- all_vars[!grepl("^(log_lik|L_rep|R_rep|L_hat|R_hat|lp__)", all_vars)]

	fn <- if (type == "median") stats::median else mean
	vals <- vapply(pars, function(p) fn(draws[[p]]), numeric(1))
	names(vals) <- pars
	vals
}

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

#' Posterior Summary for a BDEB Fit
#'
#' Returns a tidy summary table of posterior draws for model parameters
#' (and optionally derived quantities), analogous to
#' [stats::summary.lm()] for frequentist fits.
#'
#' @param object A [bdeb_fit()] object.
#' @param pars Character vector of parameter names.  Default: all model
#'   parameters (excludes `log_lik`, `L_hat`, `L_rep`, `R_hat`, `R_rep`,
#'   `lp__`, and the internal `p_Am_new`).
#' @param prob Probability for the central credible interval. Default
#'   0.90 (5th/95th percentiles).
#' @param ... Ignored.
#' @return A `posterior::draws_summary` data frame with columns
#'   `variable`, `mean`, `sd`, `median`, two quantile columns named
#'   by their percentile (e.g. `"5%"` / `"95%"` for `prob = 0.90`),
#'   `rhat`, `ess_bulk`, and `ess_tail`.
#' @export
summary.bdeb_fit <- function(object, pars = NULL, prob = 0.90, ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}

	draws <- posterior::as_draws_df(object$fit$draws())

	if (is.null(pars)) {
		all_vars <- posterior::variables(draws)
		pars <- all_vars[!grepl("^(log_lik|L_hat|L_rep|R_hat|R_rep|lp__|p_Am_new)", all_vars)]
	}

	alpha <- (1 - prob) / 2
	posterior::summarise_draws(
		posterior::subset_draws(draws, variable = pars),
		"mean", "sd", "median",
		"lower" = ~ quantile(.x, alpha),
		"upper" = ~ quantile(.x, 1 - alpha),
		"rhat",
		"ess_bulk",
		"ess_tail"
	)
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

# Internal: drop posterior metadata columns ('.chain', '.iteration',
# '.draw') from a matrix obtained via as.matrix() on a draws_df.
.bdeb_strip_meta <- function(M) {
	M[, !colnames(M) %in% c(".chain", ".iteration", ".draw"), drop = FALSE]
}

# Internal: names of model parameters (i.e. excluding log_lik, L_hat,
# L_rep, R_hat, R_rep, lp__, and the internal p_Am_new).
.bdeb_par_names <- function(draws) {
	all_vars <- posterior::variables(draws)
	all_vars[!grepl("^(log_lik|L_hat|L_rep|R_hat|R_rep|lp__|p_Am_new)", all_vars)]
}

#' Posterior Credible Intervals for BDEB Model Parameters
#'
#' Bayesian counterpart to [stats::confint()]: returns the posterior
#' central credible interval for each model parameter.
#'
#' @param object A [bdeb_fit()] object.
#' @param parm Character vector of parameter names.  Default: all model
#'   parameters.
#' @param level Probability mass of the credible interval.  Default 0.95.
#' @param ... Ignored.
#' @return A matrix with one row per parameter and two columns named by
#'   their percentile (e.g. `"2.5%"`, `"97.5%"` for `level = 0.95`).
#' @export
confint.bdeb_fit <- function(object, parm = NULL, level = 0.95, ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}
	if (length(level) != 1L || !is.numeric(level) ||
	    level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single numeric in (0, 1).")
	}

	draws <- posterior::as_draws_df(object$fit$draws())
	pars <- .bdeb_par_names(draws)
	if (!is.null(parm)) {
		missing <- setdiff(parm, pars)
		if (length(missing) > 0) {
			cli::cli_abort("Unknown parameter{?s}: {.val {missing}}.")
		}
		pars <- parm
	}

	alpha <- (1 - level) / 2
	pcts <- c(alpha, 1 - alpha)
	col_labels <- paste0(format(100 * pcts, trim = TRUE), "%")

	out <- vapply(pars, function(p) {
		stats::quantile(draws[[p]], probs = pcts, names = FALSE)
	}, numeric(2))
	out <- t(out)
	rownames(out) <- pars
	colnames(out) <- col_labels
	out
}

#' Number of Observations Used in a BDEB Fit
#'
#' Returns the total count of growth (and, where relevant, reproduction)
#' observations that contributed to the likelihood.  Matches
#' [stats::nobs()] in spirit.
#'
#' @param object A [bdeb_fit()] object.
#' @param ... Ignored.
#' @return Integer scalar: total number of observations.
#' @export
nobs.bdeb_fit <- function(object, ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}
	dat <- object$model$data
	n <- 0L
	if (!is.null(dat$growth))       n <- n + nrow(dat$growth)
	if (!is.null(dat$reproduction)) n <- n + nrow(dat$reproduction)
	as.integer(n)
}

#' Fitted Values from a BDEB Fit
#'
#' Posterior point estimate (median or mean) of the latent length
#' \eqn{\hat{L}_i} at each observation.  Bayesian counterpart of
#' [stats::fitted()].
#'
#' Currently supported for `"individual"` and `"growth_repro"` models.
#' For `"hierarchical"` and `"debtox"` use [bdeb_predict()] instead.
#'
#' @param object A [bdeb_fit()] object.
#' @param type One of `"median"` (default) or `"mean"`.
#' @param ... Ignored.
#' @return Named numeric vector of fitted values, one per observation.
#' @export
fitted.bdeb_fit <- function(object, type = c("median", "mean"), ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}
	if (!object$model$type %in% c("individual", "growth_repro")) {
		cli::cli_abort(c(
			"{.fn fitted} is not available for {.val {object$model$type}} models.",
			"i" = "Use {.fn bdeb_predict} for trajectory predictions."
		))
	}
	type <- match.arg(type)
	draws <- posterior::as_draws_df(object$fit$draws())
	all_vars <- posterior::variables(draws)
	L_hat_vars <- grep("^L_hat\\[", all_vars, value = TRUE)
	if (length(L_hat_vars) == 0L) {
		cli::cli_abort("No {.var L_hat} variables found in the fit.")
	}
	fn <- if (type == "median") stats::median else mean
	vals <- vapply(L_hat_vars, function(p) fn(draws[[p]]), numeric(1))
	names(vals) <- L_hat_vars
	vals
}

#' Residuals from a BDEB Fit
#'
#' Observed minus fitted length for each observation, using the posterior
#' point estimate from [fitted()].  Bayesian counterpart of
#' [stats::residuals()].
#'
#' Currently supported for `"individual"` and `"growth_repro"` models.
#'
#' @param object A [bdeb_fit()] object.
#' @param type Currently only `"response"` is supported (raw residuals).
#' @param ... Ignored.
#' @return Named numeric vector of residuals.
#' @export
residuals.bdeb_fit <- function(object, type = "response", ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}
	type <- match.arg(type, choices = "response")
	fits <- fitted(object)
	L_obs <- as.vector(object$model$stan_data$L_obs)
	if (length(L_obs) != length(fits)) {
		cli::cli_abort(c(
			"Length mismatch: {length(L_obs)} observations vs {length(fits)} fitted values.",
			"i" = "Cannot compute residuals for this model layout."
		))
	}
	out <- L_obs - fits
	names(out) <- names(fits)
	out
}

#' Posterior Covariance Matrix of BDEB Model Parameters
#'
#' Computes the empirical covariance matrix of the posterior draws for
#' all model parameters.  Bayesian counterpart of [stats::vcov()].
#'
#' @param object A [bdeb_fit()] object.
#' @param ... Ignored.
#' @return A symmetric numeric matrix with model parameters on rows
#'   and columns.
#' @export
vcov.bdeb_fit <- function(object, ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}
	draws <- posterior::as_draws_df(object$fit$draws())
	pars <- .bdeb_par_names(draws)
	M <- as.matrix(posterior::subset_draws(draws, variable = pars))
	M <- .bdeb_strip_meta(M)
	stats::cov(M)
}

#' Log-Likelihood of a BDEB Fit
#'
#' Computes the log-pointwise predictive density (lppd):
#' \deqn{\mathrm{lppd} = \sum_{i=1}^{N} \log\!\left(
#'   \frac{1}{S} \sum_{s=1}^{S} \exp(\log p(y_i \mid \theta^{(s)}))
#' \right),}
#' the Bayesian analogue of [stats::logLik()].  This is the natural
#' point summary of model fit reported by `bdeb_loo()` (which adds
#' Pareto-smoothed importance sampling for cross-validation).
#'
#' Requires the Stan model to store per-observation `log_lik` in
#' `generated quantities`.  Currently available for `"individual"` and
#' `"growth_repro"` models; `"hierarchical"` and `"debtox"` compute the
#' likelihood inside `reduce_sum` and do not expose individual terms.
#'
#' @param object A [bdeb_fit()] object.
#' @param ... Ignored.
#' @return A `logLik` object with attributes `df` (number of model
#'   parameters) and `nobs` (number of observations).
#' @export
logLik.bdeb_fit <- function(object, ...) {
	if (!inherits(object, "bdeb_fit")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_fit} object.")
	}
	if (object$model$type %in% c("hierarchical", "debtox")) {
		cli::cli_abort(c(
			"{.fn logLik} is not available for {.val {object$model$type}} models.",
			"i" = "Per-observation {.var log_lik} is not stored in generated quantities."
		))
	}
	draws <- posterior::as_draws_df(object$fit$draws())
	all_vars <- posterior::variables(draws)
	ll_vars <- grep("^log_lik\\[", all_vars, value = TRUE)
	if (length(ll_vars) == 0L) {
		cli::cli_abort("No {.var log_lik} variables found in the fit.")
	}

	ll <- as.matrix(posterior::subset_draws(draws, variable = ll_vars))
	ll <- .bdeb_strip_meta(ll)
	# Numerically stable lppd per observation.
	pointwise <- apply(ll, 2, function(x) {
		mx <- max(x)
		mx + log(mean(exp(x - mx)))
	})

	pars <- .bdeb_par_names(draws)
	structure(
		sum(pointwise),
		df    = length(pars),
		nobs  = length(ll_vars),
		class = "logLik"
	)
}

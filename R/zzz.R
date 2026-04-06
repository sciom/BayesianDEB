# Package load / attach hooks

.bdeb_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
	.bdeb_env$stan_dir <- system.file("stan", package = pkgname)
}

.onAttach <- function(libname, pkgname) {
	ver <- utils::packageVersion(pkgname)
	has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)

	msg <- paste0("BayesianDEB v", ver,
	              " -- Bayesian Dynamic Energy Budget Modelling")

	if (!has_cmdstanr) {
		msg <- paste0(msg, "\n",
			"Note: 'cmdstanr' not found. Model fitting requires it.\n",
			"  install.packages(\"cmdstanr\",\n",
			"    repos = c(\"https://stan-dev.r-universe.dev\",\n",
			"              getOption(\"repos\")))\n",
			"  cmdstanr::install_cmdstan()")
	}

	packageStartupMessage(msg)
}

#' Get path to a bundled Stan model file
#' @param model_name Name of the Stan model (without .stan extension).
#' @return Full path to the `.stan` file.
#' @keywords internal
stan_file <- function(model_name) {
	path <- file.path(.bdeb_env$stan_dir, paste0(model_name, ".stan"))
	if (!file.exists(path)) {
		cli::cli_abort("Stan model {.file {model_name}.stan} not found in package.")
	}
	path
}

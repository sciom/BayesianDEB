# Package load / attach hooks

.bdeb_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
	.bdeb_env$stan_dir <- system.file("stan", package = pkgname)
}

.onAttach <- function(libname, pkgname) {
	packageStartupMessage("BayesianDEB v", utils::packageVersion(pkgname))
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

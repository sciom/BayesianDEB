# =============================================================
# BayesianDEB JSS replication package — setup
# Hackenberger, Djerdj & Hackenberger (2026)
#
# Sourced at the top of every numbered script (01_, 02_, 03_).
# Defines:
#   - common libraries
#   - reproducible seed
#   - colour palette used across all manuscript figures
#   - BDEB_MODE: "lite" (default; ~30 min total) vs "full"
#     (~3 h total; matches the publication-grade settings)
#   - paths to data/ and outputs/
# =============================================================

suppressPackageStartupMessages({
	library("BayesianDEB")
	library("posterior")
	library("ggplot2")
})

# ----- Reproducibility ----------------------------------------------

set.seed(20260418L)
options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE)
options(mc.cores = max(1L, parallel::detectCores() - 1L))

# ----- Run mode -----------------------------------------------------
# BDEB_MODE = "lite" (default): chains = 2, iter = 300 + 300
#   Total runtime ~30 min on a 4-core machine.  Useful for a quick
#   sanity check or a CRAN-equivalent timing budget.
# BDEB_MODE = "full":           chains = 4, iter = 1000 + 1000
#   Reproduces the exact tables and figures of the manuscript.
#   Total runtime ~3 h on a standard 8-core machine.
# Switch mode without editing scripts:
#   Sys.setenv(BDEB_MODE = "full")
# or set BDEB_MODE=full in the calling shell.

BDEB_MODE <- Sys.getenv("BDEB_MODE", "lite")
if (!BDEB_MODE %in% c("lite", "full")) {
	stop("BDEB_MODE must be 'lite' or 'full'; got '", BDEB_MODE, "'.")
}

ITER_CFG <- if (BDEB_MODE == "full") {
	list(chains = 4L, warmup = 1000L, sampling = 1000L)
} else {
	list(chains = 2L, warmup = 300L,  sampling = 300L)
}

cat(sprintf("BDEB_MODE = %s  (chains = %d, warmup = %d, sampling = %d)\n",
            BDEB_MODE, ITER_CFG$chains, ITER_CFG$warmup, ITER_CFG$sampling))

# ----- Paths --------------------------------------------------------

REP_ROOT  <- if (interactive()) getwd() else dirname(sys.frame(1)$ofile)
if (basename(REP_ROOT) != "replication") {
	# Fallback when sourced from elsewhere
	candidate <- file.path(getwd(), "replication")
	if (dir.exists(candidate)) REP_ROOT <- candidate
}
DATA_DIR    <- file.path(REP_ROOT, "data")
OUTPUTS_DIR <- file.path(REP_ROOT, "outputs")
SBC_DIR     <- file.path(REP_ROOT, "sbc")
if (!dir.exists(OUTPUTS_DIR)) dir.create(OUTPUTS_DIR, recursive = TRUE)

# ----- Colour palette (consistent across manuscript figures) --------

col_posterior  <- "#2166AC"   # blue: posterior draws / trajectories
col_observed   <- "#B2182B"   # red:  observed data points
col_replicated <- "grey60"    # grey: PPC replicated data
col_prior      <- "#D6604D"   # light red: prior densities

# ----- Helper: timed wrapper for slow blocks ------------------------

timed <- function(label, expr) {
	t0 <- proc.time()[["elapsed"]]
	cli::cli_h2("{label}")
	res <- force(expr)
	dt <- proc.time()[["elapsed"]] - t0
	cli::cli_alert_info("{label}: {round(dt, 1)} s")
	res
}

# ----- Helper: save a ggplot to outputs/ ----------------------------

save_fig <- function(plot, name, width = 6, height = 4) {
	path <- file.path(OUTPUTS_DIR, name)
	ggsave(path, plot, width = width, height = height)
	cli::cli_alert_success("Saved {.file {path}}")
	invisible(path)
}

# ----- Sanity check -------------------------------------------------

cli::cli_alert_info("BayesianDEB version: {packageVersion('BayesianDEB')}")
cli::cli_alert_info("Output directory:    {.file {OUTPUTS_DIR}}")

# =============================================================
# BayesianDEB JSS replication — single master script
# Hackenberger, Djerdj & Hackenberger (2026)
#
# Run EVERYTHING with one command:
#
#     Rscript replicate_all.R
#
# By default this LOADS the archived posterior draws shipped in
# outputs/ and sbc/ and regenerates every figure, table and printed
# number in the manuscript.  No MCMC is run, so it completes in a few
# minutes on a laptop -- well inside the JSS one-hour budget.
#
# To refit every model from scratch instead (hours of MCMC; requires a
# working CmdStan toolchain), set:
#
#     BDEB_RECOMPUTE=true BDEB_MODE=full Rscript replicate_all.R
#
# The long simulation-based-calibration (SBC) runs are NEVER executed by
# this script; their archived rank matrices (sbc/*.rds) reproduce
# Figures 3-4 and Tables 4-5 directly.  See sbc/README.md to regenerate
# them (45 min -- 12 h per model, not part of the replication budget).
# =============================================================

t_start <- proc.time()[["elapsed"]]

.script_dir <- (function() {
	args <- commandArgs(trailingOnly = FALSE)
	m <- grep("^--file=", args, value = TRUE)
	if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
	ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
	if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
	getwd()
})()

run_section <- function(file) {
	t0 <- proc.time()[["elapsed"]]
	cat(sprintf("\n############### %s ###############\n", file))
	source(file.path(.script_dir, file), local = new.env(), chdir = TRUE)
	cat(sprintf("[%s done in %.1f s]\n", file,
	            proc.time()[["elapsed"]] - t0))
}

run_section("01_illustrations.R")   # Section 5
run_section("02_validation.R")      # Section 6 (incl. SBC figures from cache)
run_section("03_case_studies.R")    # Section 7

# Section 9 (deBInfer comparison) is a live wall-clock benchmark, so it
# only makes sense when recomputing -- caching cannot reproduce a timing.
# It runs both samplers from scratch (~10 min) and needs the deBInfer
# package.  Skipped in the default cache run.
RECOMPUTE <- tolower(Sys.getenv("BDEB_RECOMPUTE", "false")) %in%
	c("true", "1", "yes")
if (RECOMPUTE && requireNamespace("deBInfer", quietly = TRUE)) {
	run_section("comparison_debinfer.R")
} else if (!RECOMPUTE) {
	cat("\n[comparison_debinfer.R skipped in cache mode; it is a live",
	    "timing benchmark.\n set BDEB_RECOMPUTE=true to run it.]\n")
} else {
	cat("\n[comparison_debinfer.R skipped: package 'deBInfer' not installed]\n")
}

dt <- proc.time()[["elapsed"]] - t_start
cat(sprintf("\n=================================================\n"))
cat(sprintf("ALL SECTIONS COMPLETE in %.1f s (%.1f min)\n", dt, dt / 60))
cat(sprintf("Mode: %s\n",
            if (tolower(Sys.getenv("BDEB_RECOMPUTE", "false")) %in%
                c("true", "1", "yes")) "RECOMPUTE (refit from scratch)"
            else "cache (archived draws)"))
cat("Figures and tables written to replication/outputs/\n")
cat("=================================================\n")

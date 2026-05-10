# =============================================================
# Section 7: Real-data case studies
#   7.1 E. fetida growth (Neuhauser et al. 1980; gr0226 in EGrowth)
#   7.2 Cd toxicity in E. andrei (Van Gestel 1991; gr0119-gr0123)
#
# Data are bundled locally in data/curves.txt (A5 fix);
# temperature spec uses T_obs (not the legacy T = ...) -- A6 fix.
# =============================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

.script_dir <- (function() {
	args <- commandArgs(trailingOnly = FALSE)
	m <- grep("^--file=", args, value = TRUE)
	if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
	ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
	if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
	getwd()
})()
source(file.path(.script_dir, "00_setup.R"))


# -------------------------------------------------------------
# Load EGrowth curves (locally bundled — A5)
# -------------------------------------------------------------

curves_path <- file.path(DATA_DIR, "curves.txt")
if (!file.exists(curves_path)) {
	stop("Missing ", curves_path, " -- see replication/README.md for ",
	     "instructions on how to download EGrowth curves.txt locally.")
}
curves <- read.delim(curves_path)


# -------------------------------------------------------------
# 7.1 Real-data Eisenia fetida growth (Neuhauser 1980)
# -------------------------------------------------------------

cli::cli_h1("Section 7.1: Real-data Eisenia fetida (Neuhauser 1980)")

# Curve gr0226: 37 weekly observations over 250 d, 25 deg C.
neuh <- curves[curves$CURVE_ID == "gr0226", c("time", "bm")]
names(neuh) <- c("time", "wet_mass_mg")

# Convert wet mass (mg) to structural length (cm):
#   V = (W_mg / 1000) / d_V;  L = V^(1/3)
# For E. fetida: d_V ~ 1.05 g/cm^3.
d_V <- 1.05
neuh$V      <- (neuh$wet_mass_mg / 1000) / d_V
neuh$length <- neuh$V^(1 / 3)

eisenia_real <- data.frame(id = 1, time = neuh$time, length = neuh$length)

cat("\n--- Neuhauser 1980 growth data ---\n")
cat("Time points:", nrow(eisenia_real), "\n")
cat("Time range :", range(eisenia_real$time), "days\n")
cat("Length range:", round(range(eisenia_real$length), 3), "cm\n")

# Figure: raw data (fig:rawdata)
p_raw <- ggplot(eisenia_real, aes(time, length)) +
	geom_point(size = 2.5) + geom_line(alpha = 0.3) +
	theme_bw(base_size = 12) +
	labs(x = "Time (days)", y = "Structural length (cm)",
	     title = "E. fetida growth (Neuhauser et al. 1980)")
save_fig(p_raw, "fig_rawdata.pdf", width = 6, height = 4)

# Fit individual model with Arrhenius correction.
# A6 fix: use T_obs (not T = ...); the field was renamed in 0.1.4
# to avoid shadowing R's built-in T = TRUE.
dat_real <- bdeb_data(growth = eisenia_real, f_food = 1.0)
mod_real <- bdeb_model(dat_real, type = "individual",
	priors = list(
		p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
		p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
		kappa   = prior_beta(a = 3, b = 2),
		v       = prior_lognormal(mu = -1.5, sigma = 0.5),
		E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
		sigma_L = prior_halfnormal(sigma = 0.05)),
	temperature = list(T_obs = 298.15, T_ref = 293.15, T_A = 8000))

fit_real <- bdeb_fit(mod_real,
                     chains        = ITER_CFG$chains,
                     iter_warmup   = ITER_CFG$warmup,
                     iter_sampling = ITER_CFG$sampling,
                     adapt_delta   = 0.95, max_treedepth = 12,
                     seed = 42, refresh = 0)

cat("\n--- Section 7.1: Diagnostics ---\n")
diag <- bdeb_diagnose(fit_real)

# Table 6: Real data posterior
cat("\n--- Table 6: Real data posterior ---\n")
s_real <- summary(fit_real,
  pars = c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"), prob = 0.90)
print(as.data.frame(s_real), digits = 3)

# Derived quantities
cat("\n--- Derived quantities (real data) ---\n")
der_real <- bdeb_derived(fit_real,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"), f = 1.0)
print(as.data.frame(summarise_draws(der_real, "mean", "sd",
  "q5"  = ~quantile(.x, 0.05),
  "q95" = ~quantile(.x, 0.95))), digits = 3)

# Posterior correlation (Spearman)
draws_mat <- as_draws_matrix(fit_real$fit$draws())
pars_draws <- draws_mat[, c("p_Am", "p_M", "kappa", "v", "E_G",
                             "sigma_L")]
cat("\n--- Posterior correlation (Spearman) ---\n")
print(round(cor(pars_draws, method = "spearman"), 2))

# Figures: trace, posterior, pairs, trajectory, PPC, prior predictive
p_trace <- plot(fit_real, type = "trace",
                pars = c("p_Am", "p_M", "kappa", "sigma_L"))
save_fig(p_trace, "fig_trace.pdf", width = 8, height = 5)

p_post <- plot(fit_real, type = "posterior",
               pars = c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"))
save_fig(p_post, "fig_posterior.pdf", width = 8, height = 5)

if (requireNamespace("gridExtra", quietly = TRUE)) {
	p_pairs <- plot(fit_real, type = "pairs",
	                pars = c("p_Am", "p_M", "kappa", "E_G"))
	save_fig(p_pairs, "fig_pairs.pdf", width = 8, height = 8)
} else {
	cli::cli_alert_warning(
		"Package gridExtra not installed -- skipping pairs plot."
	)
}

p_traj <- plot(fit_real, type = "trajectory", n_draws = 200, seed = 1)
save_fig(p_traj, "fig_trajectory.pdf", width = 6, height = 4)

ppc <- bdeb_ppc(fit_real, type = "growth")
p_ppc <- plot(ppc, n_draws = 200)
save_fig(p_ppc, "fig_ppc.pdf", width = 6, height = 4)

pp <- bdeb_prior_predictive(mod_real, n_draws = 500, seed = 42)
p_prior <- plot(pp, n_draws = 100)
save_fig(p_prior, "fig_prior_predictive.pdf", width = 6, height = 4)

p_pp <- plot(fit_real, type = "prior_posterior")
save_fig(p_pp, "fig_prior_posterior.pdf", width = 8, height = 5)


# -------------------------------------------------------------
# 7.2 Real-data DEBtox: Cd in E. andrei (Van Gestel 1991)
# -------------------------------------------------------------
# Optional in lite mode (slow); enabled in full mode.

if (BDEB_MODE == "full") {
	cli::cli_h1("Section 7.2: Real-data DEBtox (Van Gestel 1991)")

	ids      <- c("gr0119", "gr0120", "gr0121", "gr0122", "gr0123")
	conc_map <- c(gr0119 = 0, gr0120 = 10, gr0121 = 32,
	              gr0122 = 100, gr0123 = 320)

	vg <- curves[curves$CURVE_ID %in% ids, ]
	vg$length        <- ((vg$bm / 1000) / d_V)^(1/3)
	vg$concentration <- conc_map[vg$CURVE_ID]
	vg$id            <- vg$CURVE_ID

	dat_vg <- bdeb_data(
		growth = vg[, c("id", "time", "length", "concentration")],
		concentration = conc_map,
		f_food = 1.0
	)

	mod_vg <- bdeb_tox(dat_vg, stress = "assimilation",
		priors = list(
			z_w = prior_lognormal(mu = 3.5, sigma = 1.0),
			b_w = prior_lognormal(mu = -4.0, sigma = 1.5)))

	fit_vg <- bdeb_fit(mod_vg,
	                   chains        = ITER_CFG$chains,
	                   iter_warmup   = ITER_CFG$warmup,
	                   iter_sampling = ITER_CFG$sampling,
	                   adapt_delta   = 0.95, max_treedepth = 12,
	                   threads_per_chain = 2,
	                   seed = 42, refresh = 0)

	cat("\n--- Section 7.2: Van Gestel DEBtox EC50 ---\n")
	print(bdeb_ec50(fit_vg, prob = 0.90)$summary, digits = 3)

	saveRDS(fit_vg, file.path(OUTPUTS_DIR, "fit_vangestel.rds"))
} else {
	cli::cli_alert_info(
		"Section 7.2 (Van Gestel DEBtox) skipped in lite mode; ",
		"set BDEB_MODE=full to include it."
	)
}


# -------------------------------------------------------------
# Persist Eisenia real fit + reproducibility report
# -------------------------------------------------------------

saveRDS(fit_real, file.path(OUTPUTS_DIR, "fit_real_neuhauser.rds"))

cat("\n--- Reproducibility report ---\n")
bdeb_session_info(fit_real)

writeLines(capture.output(sessionInfo()),
           file.path(OUTPUTS_DIR, "sessionInfo.txt"))

cli::cli_alert_success("Section 7 case studies: done.")

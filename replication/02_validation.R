# =============================================================
# Section 6: Validation
# Reproduces:
#   6.1 Simulation-based calibration (SBC, Talts et al. 2020)
#       — long-running auxiliary scripts in sbc/; this script
#         loads pre-computed RDS results and reports the summary.
#   6.2 Posterior correlation structure and identifiability
#       (prior-to-posterior contraction as a function of N_obs)
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
# 6.1 Simulation-based calibration (SBC)
# -------------------------------------------------------------

cli::cli_h1("Section 6.1: Simulation-based calibration")

cli::cli_alert_info("SBC runs are long ({.val 45 min} -- {.val 12 h} per model).")
cli::cli_alert_info("Auxiliary scripts in {.path sbc/} regenerate the RDS files below; this script loads them.")

sbc_paths <- list(
	individual    = file.path(SBC_DIR, "sbc_results.rds"),
	hierarchical  = file.path(SBC_DIR, "sbc_hierarchical_results.rds")
)

for (name in names(sbc_paths)) {
	p <- sbc_paths[[name]]
	if (!file.exists(p)) {
		cli::cli_alert_warning(
			"{.file {basename(p)}} not found -- run sbc/sbc_{name}.R first."
		)
		next
	}
	res <- readRDS(p)
	cat(sprintf("\n--- SBC %s ---\n", name))
	if (!is.null(res$n_total)) {
		cat(sprintf("  Replications:   %d\n", res$n_total))
		cat(sprintf("  Valid:          %d\n", res$n_valid %||% NA))
	}
	if (!is.null(res$summary)) {
		print(res$summary)
	} else {
		cli::cli_alert_info("RDS keys: {names(res)}")
	}
}


# -------------------------------------------------------------
# 6.2 Identifiability / contraction analysis
# -------------------------------------------------------------

cli::cli_h1("Section 6.2: Identifiability (Table 5 + Figure)")

# True parameters (E. fetida AmP values)
TRUE_PARS <- list(
	p_Am = 5.0, p_M = 0.5, kappa = 0.75, v = 0.2,
	E_G = 400, E0 = 1.0, L0 = 0.1, sigma_L = 0.015
)

fit_priors <- list(
	p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
	p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
	kappa   = prior_beta(a = 3, b = 2),
	v       = prior_lognormal(mu = -1.5, sigma = 0.5),
	E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
	sigma_L = prior_halfnormal(sigma = 0.05)
)

# --- Prior 90% CI widths ---
N_PRIOR <- 50000L
prior_samp <- data.frame(
	p_Am    = rlnorm(N_PRIOR, 1.5, 0.5),
	p_M     = rlnorm(N_PRIOR, -1.0, 0.5),
	kappa   = rbeta (N_PRIOR, 3, 2),
	v       = rlnorm(N_PRIOR, -1.5, 0.5),
	E_G     = rlnorm(N_PRIOR, 6.0, 0.5),
	sigma_L = abs(rnorm(N_PRIOR, 0, 0.05))
)
prior_samp$L_inf <- with(prior_samp, kappa * p_Am / p_M)
prior_samp$growth_rate <- with(prior_samp, {
	k_M <- p_M / E_G
	g   <- E_G * v / (kappa * p_Am)
	k_M * g / (3 * (1 + g))
})

all_pars <- c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L",
              "L_inf", "growth_rate")
prior_width <- vapply(prior_samp[all_pars], function(x) {
	unname(diff(quantile(x, c(0.05, 0.95))))
}, numeric(1))

cat("\n--- Prior 90% CI widths ---\n")
print(round(prior_width, 4))

# --- Simulate base trajectory ---
sim <- deb_simulate(t_max = 84, p_Am = 5, p_M = 0.5, kappa = 0.75,
                    v = 0.2, E_G = 400, E0 = 1, L0 = 0.1, f = 1,
                    dt = 0.05)

# --- Time-point configurations ---
# In lite mode use a smaller subset to keep runtime reasonable.
t_configs <- if (BDEB_MODE == "full") {
	list(
		"5"  = seq(0, 84, length.out = 5),
		"9"  = seq(0, 84, length.out = 9),
		"13" = seq(0, 84, by = 7),
		"25" = seq(0, 84, length.out = 25),
		"37" = seq(0, 84, length.out = 37)
	)
} else {
	list(
		"5"  = seq(0, 84, length.out = 5),
		"13" = seq(0, 84, by = 7),
		"37" = seq(0, 84, length.out = 37)
	)
}

# --- Fit each configuration ---
results <- list()
for (name in names(t_configs)) {
	t_obs  <- t_configs[[name]]
	L_true <- approx(sim$time, sim$L, xout = t_obs)$y

	set.seed(123)
	L_obs <- rnorm(length(t_obs), L_true, TRUE_PARS$sigma_L)

	df  <- data.frame(id = 1, time = t_obs, length = L_obs)
	dat <- bdeb_data(growth = df, f_food = 1.0)
	mod <- bdeb_model(dat, type = "individual", priors = fit_priors)

	fit <- bdeb_fit(mod,
	                chains        = ITER_CFG$chains,
	                iter_warmup   = ITER_CFG$warmup,
	                iter_sampling = ITER_CFG$sampling,
	                adapt_delta = 0.95, parallel_chains = ITER_CFG$chains,
	                refresh = 0, seed = 42)

	draws <- as.data.frame(as_draws_matrix(fit$fit$draws()))
	draws$L_inf <- draws$kappa * draws$p_Am / draws$p_M
	draws$growth_rate <- with(draws, {
		k_M <- p_M / E_G
		g   <- E_G * v / (kappa * p_Am)
		k_M * g / (3 * (1 + g))
	})

	post_width <- vapply(draws[, all_pars, drop = FALSE], function(x) {
		unname(diff(quantile(x, c(0.05, 0.95))))
	}, numeric(1))
	contraction <- pmax(1 - post_width / prior_width, 0)

	results[[name]] <- list(
		n_obs       = length(t_obs),
		post_width  = post_width,
		contraction = contraction
	)
	cli::cli_alert_success("Identifiability fit for N = {length(t_obs)}: done")
}

# --- Contraction table ---
cat("\n--- Contraction (Table 5: 1 = data-determined, 0 = prior-dominated) ---\n\n")
cat(sprintf("%-15s", "Parameter"))
for (name in names(results))
	cat(sprintf(" N=%2d", results[[name]]$n_obs))
cat("\n")
cat(strrep("-", 15 + 6 * length(results)), "\n", sep = "")
for (p in all_pars) {
	cat(sprintf("%-15s", p))
	for (name in names(results))
		cat(sprintf(" %4.2f", results[[name]]$contraction[p]))
	cat("\n")
}

# --- Contraction figure (fig:contraction) ---
cont_df <- do.call(rbind, lapply(names(results), function(name) {
	data.frame(
		n_obs       = results[[name]]$n_obs,
		parameter   = all_pars,
		contraction = results[[name]]$contraction[all_pars],
		stringsAsFactors = FALSE)
}))
par_labels <- c(
	p_Am = "p[Am]", p_M = "p[M]", kappa = "kappa",
	v = "v", E_G = "E[G]", sigma_L = "sigma[L]",
	L_inf = "L[infinity]", growth_rate = "dot(r)[B]"
)
cont_df$par_label <- factor(par_labels[cont_df$parameter],
                            levels = par_labels)
cont_df$type <- ifelse(cont_df$parameter %in% c("L_inf", "growth_rate"),
                       "Derived", "Primary")

p_cont <- ggplot(cont_df, aes(x = n_obs, y = contraction,
                              colour = par_label, linetype = type)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") +
	scale_x_continuous(breaks = sort(unique(cont_df$n_obs))) +
	scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
	scale_linetype_manual(values = c("Primary" = "solid",
	                                 "Derived" = "dashed")) +
	theme_bw(base_size = 10) +
	labs(x = "Number of observations", y = "Contraction c",
	     colour = "Parameter", linetype = "Type") +
	theme(legend.position = "right")

save_fig(p_cont, "fig_contraction.pdf", width = 7, height = 4)

saveRDS(results,
        file.path(OUTPUTS_DIR, "02_identifiability_results.rds"))

cli::cli_alert_success("Section 6 validation: done.")
